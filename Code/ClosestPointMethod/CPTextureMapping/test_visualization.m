clear all;
close all;
%% load initial mapping of image onto sphere.

load('InitialMaps/SphereMDS41.mat');
W = double(U);

%% initial planar image.
figure;
PlanarImage = reshape(W,size(xS));
imagesc(PlanarImage)

%% triangulation of surface using Delaunay triangulation.
xSg = xS(:); ySg = yS(:); zSg = zS(:);
DT = delaunayTriangulation(xSg,ySg,zSg);
Tri = freeBoundary(DT);
figure;
trisurf(Tri,xSg,ySg,zSg,W,'EdgeColor','none');

% %% scatter plot of the points with color.
% figure;
% scatter3(xSg,ySg,zSg,20,W,'fill')

