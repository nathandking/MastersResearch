%% First attempt to do manifold mapping with closest point method.
% This code does not implement the cp(x) on the plane first before the
% cp(x) onto the sphere. However since we are just doing the laplacian of
% this and the computational grid is uniform this should be the same result
% since the laplacian should not change points away from the boundary.
% However the general approach needs to be implemented.
clear all
close all
restoredefaultpath;
addpath('cp_matrices_Sphere');
tic
%%
% 3D example on a sphere
% Construct a grid in the embedding space

dx = 0.05;                   % grid size

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

%% Construct an interpolation matrix for closest point

% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(300);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

%% Create Laplacian matrix for heat equation

disp('Constructing laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);


%% load initial mapping of image onto sphere.
load('InitialMaps/test_subset_map/sphere_n/Sphere_MDS_noncp301.mat');
load('InitialMaps/test_subset_map/sphere_n/Sphere_MDS_noncp_color90601.mat');
W = double(Urgb);

%% Interpolation to get color onto computational points.

xS1 = cpX(:,1); yS1 = cpX(:,2); zS1 = cpX(:,3);
% X = [xS1, yS1, zS1];
% Y = [u1, u2, u3];
% NS = KDTreeSearcher(X);
% IDX = knnsearch(NS,Y,'K',9);
U1 = griddata(xS1, yS1, zS1, W(:,1), u1_init, u2_init, u3_init,'nearest');
U2 = griddata(xS1, yS1, zS1, W(:,2), u1_init, u2_init, u3_init,'nearest');
U3 = griddata(xS1, yS1, zS1, W(:,3), u1_init, u2_init, u3_init,'nearest');

Uc = [U1, U2, U3];
%% Visualize initial map.
% DT = delaunayTriangulation(u1(:),u2(:),u3(:));
% Tri = freeBoundary(DT);
% subplot(1,3,1);
% trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
% axis off;

figure;
%subplot(1,3,1);
scatter3(u1_init(:),u2_init(:),u3_init(:),20,Uc,'fill');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off;

%% Add noise to map.
N1 = 0.1*rand(length(u1_init),1);
N2 = 0.1*rand(length(u2_init),1);
N3 = 0.1*rand(length(u3_init),1);
u1 = u1_init + N1;
u2 = u2_init + N2;
u3 = u3_init + N3;

% map noisy data back onto sphere and assign initial u1, u2, u3.
[u1, u2, u3] = cpSphere(u1, u2, u3);

% visualize
% DT = delaunayTriangulation(u1(:),u2(:),u3(:));
% Tri = freeBoundary(DT);
% subplot(1,3,2);
% trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
%axis off;
%%
%subplot(1,3,2);
figure;
scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off;

%% Time-stepping for the heat equation

Tf = 0.005;
dt = 0.1*dx^2;
numtimesteps = ceil(Tf/dt)
error = zeros(numtimesteps,1);
% adjust for integer number of steps
dt = Tf / numtimesteps
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew1 = u1 + dt*L*u1;
    unew2 = u2 + dt*L*u2;
    unew3 = u3 + dt*L*u3;
    
    % closest point extension
%      unew1 = E*unew1;
%      unew2 = E*unew2;
%      unew3 = E*unew3;

    [u1, u2, u3] = cpSphere(unew1, unew2, unew3);
    
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

%subplot(1,3,3);
figure;
scatter3(u1,u2,u3,20,Uc,'fill');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off

%%
u1_plot = Eplot*u1;
u2_plot = Eplot*u2;
u3_plot = Eplot*u3;

figure;
scatter3(u1_plot,u2_plot,u3_plot,20,W,'fill');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off


%%
figure;
plot(1:numtimesteps,error,'k');
xlabel('$t$','Interpreter','latex','FontSize',20);
ylabel('$\|{\bf u}({\bf x},t)-{\bf u}_0({\bf x})\|_2$','Interpreter','latex','FontSize',20);

t_explicit = toc
