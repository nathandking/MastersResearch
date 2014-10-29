clear all
close all
restoredefaultpath;
addpath('cp_matrices_Sphere');
tic
%%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.1;                   % grid size

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

% find points in surface data that have zS approximately 1.
eps = 5*1e-1;
id1 = find(abs(zS1-eps) ~= 1);
id2 = find(abs(zS1+eps) ~= 1);
id = intersect(id1,id2);
not_id = setdiff(1:size(zS1,1),id); 

% find points in closest point data that have u3_init approximately 1.
cp_id1 = find(abs(u3_init-eps) ~= 1);
cp_id2 = find(abs(u3_init+eps) ~= 1);
cp_id = intersect(cp_id1,cp_id2);
not_cp_id = setdiff(1:size(u3_init,1),id);

% stereographic project of surface data and closest point data.
% V = xS1(id)./(1-zS1(id));
% W = yS1(id)./(1-zS1(id));
% 
% cpV = u1_init(cp_id)./(1 - u3_init(cp_id));
% cpW = u2_init(cp_id)./(1 - u3_init(cp_id));

V = xS1./(1-zS1);
W = yS1./(1-zS1);

cpV = u1_init./(1 - u3_init);
cpW = u2_init./(1 - u3_init);


% form matrices for nearest neighbour calculation.
X = [V, W];
Y = [cpV, cpW];

% neasest 3 neighbour calculations to create triangles with cp in it.
NS = KDTreeSearcher(X);
IDX = knnsearch(NS,Y,'K',3);

weights = zeros(3,size(IDX,1));
cp_nx = zeros(size(IDX,1),1);
cp_ny = zeros(size(IDX,1),1);
for i = 1:size(IDX,1)
    A = [V(IDX(i,:))'; W(IDX(i,:))'; [1 1 1]];
    b = [cpV(i); cpW(i); 1];
    weights(:,i) = A\b;
    cp_nx(i) = newx(IDX(i,:))'*weights(:,i);
    cp_ny(i) = newy(IDX(i,:))'*weights(:,i);
end

cp_newx = zeros(size(u1_init,1),1);
cp_newy = zeros(size(u2_init,1),1);

cp_newx(cp_id) = cp_nx;
cp_newy(cp_id) = cp_ny;

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
% figure;
% %subplot(1,3,1);
% scatter3(u1_init(:),u2_init(:),u3_init(:),20,Uc,'fill');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
% colormap('gray');
% axis off;
