# Pinecone Surface in MATLAB

![Pinecone Surface](images/pinecone_surface.png)

This project generates a pinecone-inspired surface in MATLAB by distributing repeated scale elements over a smooth base surface.

Each scale is oriented according to the local tangent frame and normal vector of the underlying geometry, with opening amplitude modulated by a curvature-based heuristic.

## Features

- MATLAB membrane-based base surface
- Local tangent-normal frame estimation
- Curvature-informed scale opening
- Custom pinecone-style brown colormap
- High-resolution PNG export

## Requirements

- MATLAB
- Image Processing Toolbox recommended for `imgaussfilt`

## Usage

Run:

```matlab
pinecone_surface
