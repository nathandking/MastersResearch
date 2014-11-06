function [newxy, XS, sub_idx] = compute_initial_map(n)
[x, y, z] = sphere(n);
xS = x(:); yS = y(:); zS = z(:);
XS = [xS, yS, zS];

%% Determine subset that image will be mapped to.
[azi, ele, ~] = cart2sph(xS, yS, zS);
azi_g = find(azi >= -pi/4); 
azi_l = find(azi <= pi/4); 
ele_g = find(ele >= -pi/4);
ele_l = find(ele <= pi/4);

azi_intersect = intersect(azi_g, azi_l);
ele_intersect = intersect(ele_g, ele_l);
sub_idx = intersect(azi_intersect, ele_intersect);

XS_subset = XS(sub_idx,:);

%xS_sub = xS(sub_idx); yS_sub = yS(sub_idx); zS_sub = zS(sub_idx);
%XS_subset = [xS_sub, yS_sub, zS_sub];

%% Determine unique closest points to do MDS on

[uni_XS, ~, ic] = unique(XS_subset,'rows'); 

%% Compute geodesic distance between pairs of points.
% This below only works for a sphere!

M = real(acos(uni_XS*uni_XS')).^2;

%% Apply multidimensional scaling to determine flattened coordinates.

uni_newxy = MDS(M);
newxy = uni_newxy(ic,:);   % store values for all points in cpX
