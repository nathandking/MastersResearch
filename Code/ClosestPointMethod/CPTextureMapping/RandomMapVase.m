%% First attempt to do manifold mapping with closest point method.
% This code does not implement the cp(x) on the plane first before the
% cp(x) onto the sphere. However since we are just doing the laplacian of
% this and the computational grid is uniform this should be the same result
% since the laplacian should not change points away from the boundary.
% However the general approach needs to be implemented.
clear all
close all
%restoredefaultpath;
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
[cpx_S, cpy_S, cpz_S, ~] = cpVase(xx,yy,zz);
% make into vectors
cpxg_S = cpx_S(:); cpyg_S = cpy_S(:); cpzg_S = cpz_S(:);

%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the torus
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx_T, cpy_T, cpz_T, dist] = cpTorus(xx,yy,zz);
% make into vectors
cpxg_T = cpx_T(:); cpyg_T = cpy_T(:); cpzg_T = cpz_T(:);


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
cpxg_S = cpxg_S(band); cpyg_S = cpyg_S(band); cpzg_S = cpzg_S(band);
cpxg_T = cpxg_T(band); cpyg_T = cpyg_T(band); cpzg_T = cpzg_T(band);
%xg = xx(band); yg = yy(band); zg = zz(band);


%% Construct a random map from Torus to Sphere.
u1_init = zeros(size(cpxg_T));
u2_init = zeros(size(cpxg_T));
u3_init = zeros(size(cpxg_T));

k = 1000;
data = [cpxg_S, cpyg_S, cpzg_S];
[U_init,idx] = datasample(data,k);

u1_init(1:k) = U_init(:,1);
u2_init(1:k) = U_init(:,2);
u3_init(1:k) = U_init(:,3);

% u1_init = rand(size(cpxg_T));
% u2_init = rand(size(cpyg_T));
% u3_init = rand(size(cpzg_T));
% 
% [u1_init, u2_init, u3_init] = cpSphere(u1_init, u2_init, u3_init);

% TH = 2*pi*rand(size(cpxg_T));
% PH = asin(-1+2*rand(size(cpxg_T)));
% [u1_init, u2_init, u3_init] = sph2cart(TH, PH, 1);

u1 = u1_init;
u2 = u2_init;
u3 = u3_init;

% visualize
% DT = delaunayTriangulation(u1(:),u2(:),u3(:));
% Tri = freeBoundary(DT);
% subplot(1,3,2);
% trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
%axis off;
Uc = zeros(size(u1));
Uc(1:k) = ones(k,1);
%subplot(1,2,1);
figure(1);
scatter3(u1,u2,u3,20,Uc,'fill');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
%colormap('jet');

%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

E_T = interp3_matrix(x1d, y1d, z1d, cpxg_T, cpyg_T, cpzg_T, p, band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation

disp('Constructing laplacian matrix');
L_T = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% Time-stepping for the heat equation

Tf = 0.5;
dt = 0.1*dx^2;
numtimesteps = ceil(Tf/dt)
error = zeros(numtimesteps,1);
% adjust for integer number of steps
dt = Tf / numtimesteps
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew1 = u1 + dt*L_T*u1;
    unew2 = u2 + dt*L_T*u2;
    unew3 = u3 + dt*L_T*u3;

    % closest point extension
    unew1 = E_T*unew1;
    unew2 = E_T*unew2;
    unew3 = E_T*unew3;
    
    [u1, u2, u3] = cpVase(unew1, unew2, unew3);
    
    figure(1); hold off
    scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
    view([1 0 0]);
    axis([-1 1 -1 1 -1 1]);
    
    
    t = kt*dt;
    error(kt) = norm([u1,u2,u3] - [u1_init, u2_init, u3_init]); 
end

%% triangulation of surface using Delaunay triangulation.
% DT = delaunayTriangulation(u1(:),u2(:),u3(:));
% Tri = freeBoundary(DT);
% subplot(1,3,3);
% trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
%axis off;

% subplot(1,2,2);
% scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
%colormap('jet');

figure(2);
plot(1:numtimesteps,error);

t_explicit = toc