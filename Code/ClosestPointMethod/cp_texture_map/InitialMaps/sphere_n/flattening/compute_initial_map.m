function [newxy, sub_idx] = compute_initial_map(n)
[x, y, z] = sphere(n);
XS = [x(:), y(:), z(:)];
clear x y z;
%% Determine subset that image will be mapped to.
[azi, ele, ~] = cart2sph(XS(:,1), XS(:,2), XS(:,3));
azi_g = find(azi >= -pi/4); 
azi_l = find(azi <= pi/4); 
ele_g = find(ele >= -pi/4);
ele_l = find(ele <= pi/4);

azi_intersect = intersect(azi_g, azi_l);
ele_intersect = intersect(ele_g, ele_l);
sub_idx = intersect(azi_intersect, ele_intersect);
clear azi_g azi_l ele_g ele_l azi_intersect ele_intersect;

XS_subset = XS(sub_idx,:);
clear XS;

%% Determine unique closest points to do MDS on
[uni_XS, ~, ic] = unique(XS_subset,'rows'); 
clear XS_subset;
%% Compute geodesic distance between pairs of points.
% This below only works for a sphere!

M = real(acos(uni_XS*uni_XS')).^2;
clear uni_XS;

%% Apply classical multidimensional scaling to flatten coordinates.
J = eye(size(M,1)) - ones(size(M,1))./(size(M,1));
B = -0.5 * J * M * J;
B = 0.5*(B+B');       % In theory B should be symmetric, so fix this.
clear J;

% Find largest eigenvalues + their eigenvectors
opts.tol = 1e-8;
[Q, L,~] = eigs(B,2,'LA',opts);
clear B;

% Extract the coordinates.
uni_newxy = [sqrt(L(1,1)).*Q(:,1), sqrt(L(2,2)).*Q(:,2)];
clear L Q;
newxy = uni_newxy(ic,:);   % store values for all points
clear uni_newxy;