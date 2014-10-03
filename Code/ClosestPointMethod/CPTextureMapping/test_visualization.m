restoredefaultpath;
addpath('cp_matrices_Sphere');
addpath('barylag2d');
%%
% 3D example on a sphere
% Construct a grid in the embedding space
tic
dx = 0.2;                   % grid size

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

%% Find unique points to apply MDS.
P = [cpxg, cpyg, cpzg];

[~, ia, ~] = unique(P,'first','rows');
ia = sort(ia);

cpxg_unique = P(ia,1);
cpyg_unique = P(ia,2);
cpzg_unique = P(ia,3);

u1 = cpxg_unique;
u2 = cpyg_unique;
u3 = cpzg_unique;

%% Apply MDS to unique list.
% Compute geodesic distance between pairs of points.
% This below only works for a sphere!
A = [u1, u2, u3];
NumPts = length(u1);
M = zeros(NumPts);
for j = 1:NumPts
    for i = 1:NumPts
        M(i,j) = (real(acos(sum(A(i,:,:).*A(j,:,:))))).^2;
    end
end

[newx, newy] = MDS(M);
toc
%%
% load('dx0p2/newx.mat');
% load('dx0p2/newy.mat');
uvColor = imread('BlurredThreeShapes512.jpg');
%uvColor = checkerboard(10,2,2);
% 
% newx = newx(band); newy = newy(band);
% newx = newx(ia); newy = newy(ia);
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(uvColor,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(uvColor,2)-1):ymax)';

[uu,vv]=meshgrid(u,v);
figure;
mesh(uu,vv,uvColor);

%find indices of the closest points in (u,v) to (newx,newy).
UV = [uu(:),vv(:)];
newXY = [newx, newy];
k = dsearchn(UV,newXY);
p = 0;
uvC = uvColor(:);
Ui = uvC(k);
%uvC = uvColor(:);
% for i = 1:length(k)
% Ui(i) = barylag2d(double(uvC(k(i)-p:k(i)+p)),...
%     uu(k(i)-p:k(i)+p),vv(k(i)-p:k(i)+p),newx(i),newy(i));
% end

figure;
imagesc(Ui)

DT = delaunayTriangulation(u1,u2,u3);
Tri = freeBoundary(DT);
figure;
trisurf(Tri,u1,u2,u3,Ui);
 figure;
 scatter3(u1,u2,u3,20,Ui)

