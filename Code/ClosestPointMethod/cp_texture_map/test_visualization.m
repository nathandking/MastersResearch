clear all;
close all;

% mimic computational grid.
dx = 0.2;
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);

[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

%% load initial mapping of image onto sphere.

load('InitialMaps/tests/old_initial_maps/SphereMDS41.mat');
W = double(U);

%% initial planar image.
figure;
PlanarImage = reshape(W,size(xS));
imagesc(PlanarImage)

%% triangulation of surface using Delaunay triangulation.
 xSg = xS(:); ySg = yS(:); zSg = zS(:);
% DT = delaunayTriangulation(xSg,ySg,zSg);
% Tri = freeBoundary(DT);
% figure;
% trisurf(Tri,xSg,ySg,zSg,W,'EdgeColor','none');

% %% scatter plot of the points with color.
% figure;
% scatter3(xSg,ySg,zSg,20,W,'fill')

X = [xSg; xx(:)];
Y = [ySg; yy(:)];
Z = [zSg; zz(:)];
colors = [W; zeros(size(xx(:),1),1)];
Vol = sqrt(X.^2+Y.^2+Z.^2)-1;
figure;
[F, V, C] = isosurface(X,Y,Z,Vol,0,colors);
patch('Vertices', V, 'Faces', F, ... 
    'FaceVertexCData', C, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');



