function pinecone_surface(nU,nV,scaleLength,scaleWidth,exportImage)
% PINECONE_SURFACE Generate a pinecone-inspired geometric surface in MATLAB.
%
% This function creates a smooth base surface and distributes pinecone-like
% scale elements over it. Each scale is oriented according to the local
% tangent frame and surface normal of the underlying geometry.
%
% Inputs (all optional):
%   nU           - number of scale columns
%   nV           - number of scale rows
%   scaleLength  - scale length
%   scaleWidth   - scale width
%   exportImage  - true/false to export PNG image
%
% Example:
%   pinecone_surface
%   pinecone_surface(22,18,0.18,0.10,true)
%
% Author: Gianluigi Riccardi
% Year: 2026

    if nargin < 1 || isempty(nU), nU = 22; end
    if nargin < 2 || isempty(nV), nV = 18; end
    if nargin < 3 || isempty(scaleLength), scaleLength = 0.18; end
    if nargin < 4 || isempty(scaleWidth), scaleWidth = 0.10; end
    if nargin < 5 || isempty(exportImage), exportImage = true; end

    close all; clc;

    % ------------------------------------------------------------
    % BASE SURFACE GRID
    % ------------------------------------------------------------
    NuSurf = 220;
    NvSurf = 180;

    u = linspace(-1.2, 1.2, NuSurf);
    v = linspace(-1.0, 1.0, NvSurf);
    [U,V] = meshgrid(u,v);

    % Smooth membrane-like base surface
    Z = 0.45*exp(-0.9*(U.^2 + 1.3*V.^2)) ...
      + 0.10*sin(2.4*pi*U).*cos(1.6*pi*V) ...
      - 0.06*(U.^2) ...
      - 0.03*(V.^2);

    % Slight smoothing for visual cleanliness
    Z = localSmooth2D(Z);

    % ------------------------------------------------------------
    % LOCAL DIFFERENTIAL GEOMETRY
    % ------------------------------------------------------------
    du = u(2) - u(1);
    dv = v(2) - v(1);

    [Zu,Zv] = gradient(Z,du,dv);

    % Normal field N = normalize([-dZ/du, -dZ/dv, 1])
    Nx = -Zu;
    Ny = -Zv;
    Nz = ones(size(Z));
    Nnorm = sqrt(Nx.^2 + Ny.^2 + Nz.^2);
    Nx = Nx ./ Nnorm;
    Ny = Ny ./ Nnorm;
    Nz = Nz ./ Nnorm;

    % Tangent directions
    Tu = cat(3, ones(size(Z)), zeros(size(Z)), Zu);
    Tv = cat(3, zeros(size(Z)), ones(size(Z)), Zv);

    Tu = normalizeField(Tu);
    Tv = normalizeField(Tv);

    % Curvature-like heuristic from Laplacian magnitude
    [Zuu,~] = gradient(Zu,du,dv);
    [~,Zvv] = gradient(Zv,du,dv);
    curvatureProxy = abs(Zuu + Zvv);
    curvatureProxy = curvatureProxy ./ max(curvatureProxy(:) + eps);

    % ------------------------------------------------------------
    % FIGURE SETUP
    % ------------------------------------------------------------
    fig = figure('Color','w','Name','Pinecone Surface in MATLAB');
    ax = axes('Parent',fig);
    hold(ax,'on');
    axis(ax,'equal');
    axis(ax,'off');
    view(ax,34,28);
    camlight(ax,'headlight');
    camlight(ax,'right');
    lighting(ax,'gouraud');
    material(ax,[0.45 0.7 0.25]);

    % Base surface
    surf(ax,U,V,Z, ...
        'EdgeColor','none', ...
        'FaceColor',[0.77 0.63 0.43], ...
        'FaceAlpha',0.92);

    % ------------------------------------------------------------
    % SCALE DISTRIBUTION
    % ------------------------------------------------------------
    uCenters = linspace(-0.95,0.95,nU);
    vCenters = linspace(-0.82,0.82,nV);

    % Alternate offset for more organic arrangement
    for j = 1:nV
        for i = 1:nU

            uu = uCenters(i);
            vv = vCenters(j);

            if mod(j,2) == 0
                uu = uu + 0.5*(uCenters(min(2,end)) - uCenters(1));
            end

            if uu < min(u) || uu > max(u) || vv < min(v) || vv > max(v)
                continue;
            end

            % Interpolate position
            zc = interp2(U,V,Z,uu,vv,'linear',NaN);
            if isnan(zc)
                continue;
            end

            % Interpolate normal and tangents
            n = [
                interp2(U,V,Nx,uu,vv,'linear',NaN)
                interp2(U,V,Ny,uu,vv,'linear',NaN)
                interp2(U,V,Nz,uu,vv,'linear',NaN)
            ];

            tu = [
                interp2(U,V,Tu(:,:,1),uu,vv,'linear',NaN)
                interp2(U,V,Tu(:,:,2),uu,vv,'linear',NaN)
                interp2(U,V,Tu(:,:,3),uu,vv,'linear',NaN)
            ];

            tv = [
                interp2(U,V,Tv(:,:,1),uu,vv,'linear',NaN)
                interp2(U,V,Tv(:,:,2),uu,vv,'linear',NaN)
                interp2(U,V,Tv(:,:,3),uu,vv,'linear',NaN)
            ];

            cval = interp2(U,V,curvatureProxy,uu,vv,'linear',NaN);

            if any(isnan([n;tu;tv;cval]))
                continue;
            end

            n  = normalizeVec(n);
            tu = normalizeVec(tu);
            tv = normalizeVec(tv);

            % Make tangents orthogonal to normal
            tu = normalizeVec(tu - dot(tu,n)*n);
            tv = normalizeVec(cross(n,tu));

            % Opening angle increases with curvature heuristic
            openAngle = deg2rad(10 + 38*cval);

            % Slight size modulation
            localLength = scaleLength * (0.90 + 0.25*cval);
            localWidth  = scaleWidth  * (0.92 + 0.18*(1-cval));

            % Organic azimuth bias
            azBias = 0.12*sin(3.1*uu + 1.7*vv);
            tuRot = cos(azBias)*tu + sin(azBias)*tv;
            tvRot = -sin(azBias)*tu + cos(azBias)*tv;

            % Create one scale patch
            [Xs,Ys,Zs] = makeScalePatch( ...
                [uu;vv;zc], ...
                tuRot, tvRot, n, ...
                localLength, localWidth, openAngle);

            % Brownish tone variation
            baseColor = [0.43 0.28 0.14];
            lightGain = 0.15 + 0.30*cval;
            faceColor = min(baseColor + lightGain*[0.55 0.40 0.18], 1);

            surf(ax,Xs,Ys,Zs, ...
                'EdgeColor','none', ...
                'FaceColor',faceColor, ...
                'FaceLighting','gouraud', ...
                'AmbientStrength',0.42, ...
                'DiffuseStrength',0.75, ...
                'SpecularStrength',0.08);
        end
    end

    % ------------------------------------------------------------
    % FINAL LOOK
    % ------------------------------------------------------------
    title(ax,'Pinecone-Inspired Surface in MATLAB', ...
        'FontWeight','bold', ...
        'FontSize',16, ...
        'Color',[0.18 0.12 0.08]);

    colormap(ax,brownMap(256));

    xlim(ax,[min(u) max(u)]);
    ylim(ax,[min(v) max(v)]);
    zlim(ax,[min(Z(:))-0.1, max(Z(:))+0.45]);

    set(fig,'Renderer','opengl');

    % ------------------------------------------------------------
    % EXPORT
    % ------------------------------------------------------------
    if exportImage
        drawnow;
        exportgraphics(fig,'pinecone_surface_example.png','Resolution',300);
        fprintf('Image exported: pinecone_surface_example.png\n');
    end
end

% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function A = localSmooth2D(A)
    K = [1 2 1; 2 4 2; 1 2 1] / 16;
    A = conv2(A,K,'same');
end

function F = normalizeField(F)
    nrm = sqrt(F(:,:,1).^2 + F(:,:,2).^2 + F(:,:,3).^2) + eps;
    F(:,:,1) = F(:,:,1) ./ nrm;
    F(:,:,2) = F(:,:,2) ./ nrm;
    F(:,:,3) = F(:,:,3) ./ nrm;
end

function v = normalizeVec(v)
    v = v(:);
    n = norm(v);
    if n < eps
        v = [0;0;1];
    else
        v = v / n;
    end
end

function [X,Y,Z] = makeScalePatch(center,tu,tv,n,L,W,openAngle)
    % Local scale coordinates
    ns = 18;
    nt = 10;

    s = linspace(0,1,ns);
    t = linspace(-1,1,nt);
    [S,T] = meshgrid(s,t);

    % Planform: narrow tip, wider root
    widthProfile = W * (0.12 + 0.88*(1 - S).^0.85);
    Xl = T .* widthProfile;
    Yl = -L * S;

    % Curvature / opening shape
    arch = 0.10*L*(1 - (T).^2).*(1 - 0.2*S);
    lift = (sin(openAngle)) * (0.55*L) .* (S.^1.4) .* (1 - 0.15*T.^2);
    Zl = arch + lift;

    % Slight longitudinal camber
    Yl = Yl + 0.08*L*(T.^2).*(S.^1.1);

    % Local basis: tv = width direction, tu = length direction
    P = center(:)' ...
      + Xl(:)*tv(:)' ...
      + Yl(:)*tu(:)' ...
      + Zl(:)*n(:)';

    X = reshape(P(:,1),size(Xl));
    Y = reshape(P(:,2),size(Xl));
    Z = reshape(P(:,3),size(Xl));
end

function cmap = brownMap(m)
    if nargin < 1
        m = 256;
    end

    anchors = [
        0.14 0.08 0.04
        0.24 0.15 0.07
        0.36 0.23 0.11
        0.52 0.35 0.17
        0.66 0.50 0.27
        0.80 0.67 0.44
    ];

    x = linspace(0,1,size(anchors,1));
    xi = linspace(0,1,m);

    cmap = zeros(m,3);
    for k = 1:3
        cmap(:,k) = interp1(x,anchors(:,k),xi,'pchip');
    end
    cmap = max(0,min(1,cmap));
end
