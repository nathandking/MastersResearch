clear all
close all
restoredefaultpath;
addpath('cp_matrices_Sphere');
tic
%%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.025;                   % grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);

%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);


%% Banding: do calculation in a narrow band around the sphere
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
xg = xx(band); yg = yy(band); zg = zz(band);


%% Function u in the embedding space
% this makes u into a vector, containing only points in the band
u1_init = cpxg;
u2_init = cpyg;
u3_init = cpzg;

%% load initial mapping of image onto sphere.

load('InitialMaps/StJohnsSphereMDS101.mat');
W = double(U);

%% Interpolation to get color onto computational points.
xS1 = xS(:); yS1 = yS(:); zS1 = zS(:);

% form matrices for nearest neighbour calculation.
X = [xS1, yS1, zS1];
[X,ia,ic] = unique(X,'rows');
Y = [u1_init, u2_init, u3_init];

% neasest 4 neighbour calculations to create triangles with cp in it.
NS = KDTreeSearcher(X);
[IDX, Dist] = knnsearch(NS,Y,'K',4);

not_near_idx = (Dist(:,1)>=1e-3);

% project 4 nearest neighbours into tangent plane of Y. 
[proj_S1, ~] = tan_plane_proj_sphere(X(IDX(:,1),:),Y);
[proj_S2, ~] = tan_plane_proj_sphere(X(IDX(:,2),:),Y);
[proj_S3, ~] = tan_plane_proj_sphere(X(IDX(:,3),:),Y);
[proj_S4, Y_3d] = tan_plane_proj_sphere(X(IDX(:,4),:),Y);

% change basis so that 3d plane is a 2d plane.
%[proj_S1_2d, proj_S2_2d, proj_S3_2d, proj_S4_2d, Y_2d] = ...
 %   plane_3d_to_2d(proj_S1, proj_S2, proj_S3, proj_S4, Y_3d);
 
proj_S1_2d = proj_S1(:,1:2);
proj_S2_2d = proj_S2(:,1:2);
proj_S3_2d = proj_S3(:,1:2);
proj_S4_2d = proj_S4(:,1:2);
Y_2d = Y_3d(:,1:2);

in = zeros(size(Y_2d,1),1);
for i = 1:size(Y_2d,1)
in(i) = inpolygon(Y_2d(i,1),Y_2d(i,2),[proj_S1_2d(i,1), proj_S2_2d(i,1), ...
    proj_S3_2d(i,1), proj_S4_2d(i,1)],[proj_S1_2d(i,1), proj_S2_2d(i,1), ...
    proj_S3_2d(i,1), proj_S4_2d(i,1)]);
end

newxy = [newx(ia),newy(ia)];
% corresponding newx and newy values are.
S1_newxy = newxy(IDX(:,1),:);
S2_newxy = newxy(IDX(:,2),:);
S3_newxy = newxy(IDX(:,3),:);
S4_newxy = newxy(IDX(:,4),:);

c_idx_123 = test_not_collinearity(proj_S1_2d, proj_S2_2d, proj_S3_2d);
c_idx_124 = test_not_collinearity(proj_S1_2d, proj_S2_2d, proj_S4_2d);
c_idx_134 = test_not_collinearity(proj_S1_2d, proj_S3_2d, proj_S4_2d);
c_idx_234 = test_not_collinearity(proj_S2_2d, proj_S3_2d, proj_S4_2d);

c_idx_125 = test_not_collinearity(proj_S1_2d, proj_S2_2d, Y_2d);
c_idx_135 = test_not_collinearity(proj_S1_2d, proj_S3_2d, Y_2d);
c_idx_145 = test_not_collinearity(proj_S1_2d, proj_S4_2d, Y_2d);
c_idx_235 = test_not_collinearity(proj_S2_2d, proj_S3_2d, Y_2d);
c_idx_245 = test_not_collinearity(proj_S2_2d, proj_S4_2d, Y_2d);
c_idx_345 = test_not_collinearity(proj_S3_2d, proj_S4_2d, Y_2d);





not_cc_idx = c_idx_123 + c_idx_124 + c_idx_134 + c_idx_234 + c_idx_125 + c_idx_135 + c_idx_145 + c_idx_235 + c_idx_245 + c_idx_345;
not_c_idx = (not_cc_idx < 10);
not_idx = find(logical(not_near_idx) +not_c_idx + logical(in) == 3);

[Y_newxy_temp] = projective_mapping(Y_2d(not_idx,:), proj_S1_2d(not_idx,:), proj_S2_2d(not_idx,:),...
    proj_S3_2d(not_idx,:),proj_S4_2d(not_idx,:), S1_newxy(not_idx,:), S2_newxy(not_idx,:), S3_newxy(not_idx,:), S4_newxy(not_idx,:));

Y_newxy = S1_newxy;
Y_newxy(not_idx,:) = Y_newxy_temp;
YNan = find(isnan(Y_newxy(:,1)));
Y_newxy(YNan,:) = S1_newxy(YNan,:);


T = imread('StJohns.jpg');
T = rgb2gray(T);
Tg = double(T(:));


% scale texture image to lie within (newx,newy) coordinates.
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

[uu,vv] = ndgrid(u,v);

%find indices of the closest points in (u,v) to (newx,newy).
% uv = [uu(:),vv(:)];
% newxy = [newx, newy];
% k = dsearchn(uv,newxy);
Uc = griddata(uu(:), vv(:), Tg, Y_newxy(:,1), Y_newxy(:,2),'natural');

% Uc = griddata(xS1, yS1, zS1, W, u1_init, u2_init, u3_init,'nearest');
% 
% %% Visualize initial map.
% % DT = delaunayTriangulation(u1(:),u2(:),u3(:));
% % Tri = freeBoundary(DT);
% % subplot(1,3,1);
% % trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
% % view([0 1 0]);
% % axis([-1 1 -1 1 -1 1]);
% % axis off;
% 
figure;
%subplot(1,3,1);
scatter3(u1_init(:),u2_init(:),u3_init(:),20,Uc,'fill');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off;
toc
