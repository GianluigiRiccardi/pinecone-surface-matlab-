%% MathWorks_Pinecone_v2.m
% Pinecone-Inspired Surface on a Membrane Geometry
%
% This script generates a pinecone-like surface by distributing oriented
% scale elements over a smooth base surface. Each scale follows the local
% curvature of the underlying geometry and is oriented according to the
% local tangent frame.
%
% The goal of the experiment is to explore nature-inspired surface tiling
% using MATLAB geometry processing and visualization tools.
%
% Author: Gianluigi Riccardi
% Year: 2026
% -------------------------------------------------------------------------
% MAIN PARAMETERS (user-tunable)
% -------------------------------------------------------------------------

close all; clc;

nU       = 22;      % number of scales along U direction
nV       = 28;      % number of scales along V direction
vps      = 13;      % vertices per scale (lateral resolution)
vpr      = 20;      % vertices per scale (depth resolution)

lift_max = 0.22;    % maximum scale opening
lift_min = 0.06;    % minimum scale opening

scale_w  = 0.055;   % half width of scale in tangent plane
scale_d  = 0.072;   % root-to-tip scale depth
tilt_deg = 28;      % scale tilt angle (degrees)

% -------------------------------------------------------------------------
% 1. BASE SURFACE (MathWorks membrane example)
% -------------------------------------------------------------------------

try
    Zraw = membrane(1,48);
catch
    % fallback surface if membrane is unavailable
    Zraw = peaks(97);
end

Ng   = size(Zraw,1);

ulin = linspace(-1,1,Ng);
vlin = linspace(-1,1,Ng);

[Ug,Vg] = meshgrid(ulin,vlin); %#ok<ASGLU>

% normalize surface height
Zraw = Zraw / max(abs(Zraw(:)));

% smooth surface slightly
Zg = imgaussfilt(Zraw,1.2);

% -------------------------------------------------------------------------
% 2. SURFACE GRADIENT AND NORMAL ESTIMATION
% -------------------------------------------------------------------------

[dZdu,dZdv] = gradient(Zg,ulin,vlin);

F_Z    = griddedInterpolant({vlin,ulin},Zg,'cubic','none');
F_dZdu = griddedInterpolant({vlin,ulin},dZdu,'cubic','none');
F_dZdv = griddedInterpolant({vlin,ulin},dZdv,'cubic','none');

% -------------------------------------------------------------------------
% 3. SCALE ROOT GRID
% -------------------------------------------------------------------------

u_roots = linspace(-0.92,0.92,nU);
v_roots = linspace(-0.92,0.92,nV);

[Ur,Vr] = meshgrid(u_roots,v_roots);

Ur = Ur(:);
Vr = Vr(:);

Zr    = F_Z(Vr,Ur);
dZr_u = F_dZdu(Vr,Ur);
dZr_v = F_dZdv(Vr,Ur);

valid = ~isnan(Zr);

Ur      = Ur(valid);
Vr      = Vr(valid);
Zr      = Zr(valid);
dZr_u   = dZr_u(valid);
dZr_v   = dZr_v(valid);

nscales = numel(Ur);

% -------------------------------------------------------------------------
% 4. LOCAL TANGENT FRAME FOR EACH SCALE
% -------------------------------------------------------------------------

Tu = [ones(nscales,1) zeros(nscales,1) dZr_u];
Tv = [zeros(nscales,1) ones(nscales,1) dZr_v];

Tu = Tu ./ (sqrt(sum(Tu.^2,2)) + 1e-9);
Tv = Tv ./ (sqrt(sum(Tv.^2,2)) + 1e-9);

Nrm = cross(Tu,Tv,2);
Nrm = Nrm ./ (sqrt(sum(Nrm.^2,2)) + 1e-9);

% ensure upward normal orientation
flip_idx = Nrm(:,3) < 0;
Nrm(flip_idx,:) = -Nrm(flip_idx,:);

% recompute orthogonal tangent
Tv = cross(Nrm,Tu,2);
Tv = Tv ./ (sqrt(sum(Tv.^2,2)) + 1e-9);

% -------------------------------------------------------------------------
% 5. LOCAL CURVATURE ESTIMATION
% -------------------------------------------------------------------------

[d2Zdu2,~]  = gradient(dZdu,ulin,vlin);
[~,d2Zdv2]  = gradient(dZdv,ulin,vlin);

F_kuu = griddedInterpolant({vlin,ulin},d2Zdu2,'linear','none');
F_kvv = griddedInterpolant({vlin,ulin},d2Zdv2,'linear','none');

kuu = F_kuu(Vr,Ur);
kvv = F_kvv(Vr,Ur);

curv_mag  = min(abs(kuu) + abs(kvv),5);
curv_norm = curv_mag / (max(curv_mag) + 1e-6);

% curvature-based scale opening
lift_val = lift_min + (lift_max - lift_min) * (1 - curv_norm);

% -------------------------------------------------------------------------
% 6. SCALE GEOMETRY CONSTRUCTION
% -------------------------------------------------------------------------

s_lat = linspace(-1,1,vps);
s_dep = linspace(0,1,vpr)';

lat_profile = 1 - s_lat.^2;

tilt_rad = deg2rad(tilt_deg);

all_X = zeros(vpr,vps*nscales);
all_Y = zeros(vpr,vps*nscales);
all_Z = zeros(vpr,vps*nscales);
all_C = zeros(vpr,vps*nscales);

for k = 1:nscales

    cidx = (k-1)*vps + (1:vps);

    tu = Tu(k,:);
    tv = Tv(k,:);
    n  = Nrm(k,:);

    sl = lift_val(k);

    open_dir = n*cos(tilt_rad) + tv*sin(tilt_rad);
    open_dir = open_dir/(norm(open_dir)+1e-9);

    lat_world  = scale_w .* s_lat;
    dep_world  = scale_d .* s_dep;
    lift_world = sl .* (s_dep.^1.5) .* lat_profile;

    Xk = Ur(k) + tu(1).*lat_world + tv(1).*dep_world + open_dir(1).*lift_world;
    Yk = Vr(k) + tu(2).*lat_world + tv(2).*dep_world + open_dir(2).*lift_world;
    Zk = Zr(k) + tu(3).*lat_world + tv(3).*dep_world + open_dir(3).*lift_world;

    all_X(:,cidx) = Xk;
    all_Y(:,cidx) = Yk;
    all_Z(:,cidx) = Zk;

    height_c = (Zr(k)+1)/2;
    depth_c  = s_dep .* ones(1,vps);

    all_C(:,cidx) = 0.4*depth_c + 0.6*height_c;

end

% -------------------------------------------------------------------------
% 7. RENDERING
% -------------------------------------------------------------------------

fig = figure('Color',[0.05 0.02 0], ...
             'Units','normalized', ...
             'Position',[0.08 0.05 0.55 0.85], ...
             'Name','MathWorks Pinecone v2');

surface(all_X,all_Y,all_Z,all_C,'EdgeColor','none');
shading interp

% pinecone-style brown colormap
stops_hex = ["#1A0D00"
             "#3B1F00"
             "#5A2E0C"
             "#7A3F14"
             "#9A5526"
             "#C48A55"];

stops_rgb = validatecolor(stops_hex,'multiple');
cmap = interp1(linspace(1,256,6),stops_rgb,1:256);

colormap(cmap)
clim([0.05 0.95])

material([0.35 0.85 0.55 25 0.25])

lighting gouraud

light('Position',[ 1.8  1.2  2.5],'Style','infinite','Color',[1.00 0.88 0.60])
light('Position',[-2.0 -0.8  1.5],'Style','infinite','Color',[0.20 0.18 0.30])
light('Position',[ 0.3  0.5 -1.8],'Style','infinite','Color',[0.10 0.08 0.04])

% -------------------------------------------------------------------------
% 8. VIEW SETTINGS
% -------------------------------------------------------------------------

ax = gca;

set(ax,'DataAspectRatio',[1 1 1], ...
       'Color',[0.05 0.02 0], ...
       'Position',[0.02 0.02 0.96 0.94], ...
       'Clipping','off')

axis off
axis([-1.15 1.15 -1.15 1.15 -1.3 1.5])

view(-37.5,30)
camzoom(1.15)

% -------------------------------------------------------------------------
% 9. EXPORT IMAGE
% -------------------------------------------------------------------------

out_path = fullfile(fileparts(mfilename('fullpath')),'MathWorks_Pinecone_v2.png');

exportgraphics(fig,out_path, ...
               'Resolution',300, ...
               'BackgroundColor',[0.05 0.02 0])

fprintf('Export completed: %s\n',out_path);
