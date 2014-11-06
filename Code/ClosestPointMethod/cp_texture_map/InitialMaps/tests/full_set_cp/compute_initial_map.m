function [newxy, cpX] = compute_initial_map(dx)
%% Construct a grid in the embedding space

x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

%% Find closest points on the surface.

[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

%% Banding: do calculation in a narrow band around the sphere.

dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order

bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

%% Store closest points in the band.

cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
cpX = [cpxg, cpyg, cpzg];

%% Determine unique closest points to do MDS on

[uni_cpX, ~, ic] = unique(cpX,'rows'); 

%% Compute geodesic distance between pairs of points.
% This below only works for a sphere!

NumPts = size(uni_cpX,1);
M = zeros(NumPts);
for j = 1:NumPts
    for i = 1:NumPts
        M(i,j) = (real(acos(sum(uni_cpX(i,:,:).*uni_cpX(j,:,:))))).^2;
    end
end

%% Apply multidimensional scaling to determine flattened coordinates.

uni_newxy = MDS(M);
newxy = uni_newxy(ic,:);   % store values for all points in cpX
