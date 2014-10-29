clear all
close all
restoredefaultpath;
addpath('../cp_matrices_Sphere');
tic
%% 
% 3D example on a sphere
% Construct a computational grid in the embedding space

dx = 0.025;                   % grid size

% make vectors of x, y, z, positions of the computational grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

nx = length(x1d);
ny = length(y1d);
nz = length(z1d);

%% Find closest points on the surface
% For each point (x,y,z), store closest point on the sphere (cpx,cpy,cpz)

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

%% load initial mapping of image onto sphere.
% gives (xS, yS, zS) corresponding (newx, newy) and color U.
load('InitialMaps/StJohnsSphereMDS101.mat');
W = double(U);

%% Interpolation to get color onto computational points. Projective mapping

%% form matrices for nearest neighbour calculation.
S = [xS(:), yS(:), zS(:)];
[XS,ia,ic] = unique(S,'rows');   % compute unique set of points
s = newx(ia); t = newy(ia);
st = [s, t];
cpX = [cpxg, cpyg, cpzg];
n_cpX = size(cpX,1);
cpX_2d = zeros(n_cpX,2);

%% find all points (xS, yS, zS) within distance r of (cpx, cpy, cpz).
r = 1e-1;    
[idx,D] = rangesearch(XS,cpX,r);

% HAVE TO CHECK TO MAKE SURE r GIVES AT LEAST 4 POINTS.
for i = 1:n_cpX
    if (size(idx{i},2)<=4)
        error('Less than four points within distance r of query point');
    end
end

%% check if closest (xS,yS,zS) are on top of (cpx, cpy, cpz).
Dist = zeros(n_cpX,1);
for i = 1:n_cpX
    Dist(i) = D{i}(1);
end

near_dist = 1e-6;
near_idx = idx(Dist <= near_dist);      % points that are within near_dist

for i = 1:size(near_idx,1)
cpX_2d(near_idx{i}(1),:) = [s(near_idx{i}(1)), t(near_idx{i}(1))];
end

cp_not_near_idx = find(Dist > near_dist);
not_near_idx = idx(Dist > near_dist);   % points that are not within near_dist

%% create non collinear and enclosing patches for each cpX.
patches = cell(size(not_near_idx,1),1);
cpX_tan_2d = cell(size(not_near_idx,1),1);
st_patch = cell(size(not_near_idx,1),1);
for i = 1:size(not_near_idx,1)
    [patches{i}, st_patch{i}, cpX_tan_2d{i}] = non_collinear_enclosing_patches(XS(not_near_idx{i},:),cpX(cp_not_near_idx(i),:),st(not_near_idx{i},:));
end

cpX_st_temp = cell(size(not_near_idx,1),1);
for i = 1:size(not_near_idx,1)
    if size(patches{i},2)==2
        cpX_st_temp{i} = projective_mapping2(cpX_tan_2d{i}, patches{i}, st_patch{i});
    end
end    

cpX_st = cell2mat(cpX_st_temp);
% need index to original closest points, just now since we do not have all
% perfect patches, so we are leaving some points out. Just needed to
% visualize the points that have good patches.

[nrows, ncols] = cellfun(@size, patches);
idxx = cp_not_near_idx(ncols == 2);


%%
% Y_newxy = S1_newxy;
% Y_newxy(not_idx,:) = Y_newxy_temp;
% YNan = find(isnan(Y_newxy(:,1)));
% Y_newxy(YNan,:) = S1_newxy(YNan,:);
% 

%%
T = imread('../Images/StJohns.jpg');
T = rgb2gray(T);
Tg = double(T(:));


% scale texture image to lie within (newx,newy) coordinates.
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

[uu,vv] = ndgrid(u,v);

Uc = griddata(uu(:), vv(:), Tg, cpX_st(:,1), cpX_st(:,2),'natural');

figure;
%subplot(1,3,1);
scatter3(cpxg(idxx),cpyg(idxx),cpzg(idxx),20,Uc,'fill');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off;
 toc
