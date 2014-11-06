function [newxy, cpX, sub_idx] = compute_initial_map(dx)
%% Construct a grid in the embedding space

x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

%% Find closest points on the surface.

[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);

cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

% figure;
% scatter3(cpxg,cpyg,cpzg,10,'fill');
% 
% axis([-1 1 -1 1 -1 1]);
% colormap('gray');
% 
% 
% newxy = []; cpX = [];
%% Banding: do calculation in a narrow band around the sphere.

dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order

bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

%% Store closest points in the band.

cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
cpX = [cpxg, cpyg, cpzg];

%% Determine subset that image will be mapped to.
[azi, ele, ~] = cart2sph(cpxg, cpyg, cpzg);
azi_g = find(azi >= -pi/4); 
azi_l = find(azi <= pi/4); 
ele_g = find(ele >= -pi/4);
ele_l = find(ele <= pi/4);

azi_intersect = intersect(azi_g, azi_l);
ele_intersect = intersect(ele_g, ele_l);
sub_idx = intersect(azi_intersect, ele_intersect);

cpxg = cpxg(sub_idx); cpyg = cpyg(sub_idx); cpzg = cpzg(sub_idx);

cpX_subset = [cpxg, cpyg, cpzg];

%% Determine unique closest points to do MDS on

[uni_cpX, ~, ic] = unique(cpX_subset,'rows'); 

%% Compute geodesic distance between pairs of points.
% This below only works for a sphere!

M = real(acos(uni_cpX*uni_cpX')).^2;

%% Apply multidimensional scaling to determine flattened coordinates.

uni_newxy = MDS(M);
newxy = uni_newxy(ic,:);   % store values for all points in cpX
