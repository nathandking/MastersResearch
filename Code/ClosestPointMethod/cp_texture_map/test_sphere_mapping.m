%% Harmonic mapping from a plane to sphere.
% visualization of the mapping is shown by diffusing a noisy texture mapped
% image. The initial texture map is computed using multidimensional
% scaling and is implemented in InitialMaps/sphere_n/.
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
[xp,yp,zp] = sphere(300);   % number of points must match initial map size.
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

%% Create Laplacian matrix for heat equation

disp('Constructing laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);


%% load initial mapping of image onto sphere.
load('InitialMaps/sphere_n/Sphere_MDS_noncp301.mat');
load('InitialMaps/sphere_n/mandrill_Sphere_MDS_noncp_color301.mat');
W = double(Urgb);

%% plot the initial texture mapped image.
U_plot_init = [Eplot*u1_init, Eplot*u2_init, Eplot*u3_init];
[V, F, C] = tri_color(U_plot_init, W); 

figure(1);
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'interp','edgecolor', 'none');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
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

%% plot the noisy texture mapped initial image.
U_plot_noisy = [Eplot*u1, Eplot*u2, Eplot*u3];
[V, F, C] = tri_color(U_plot_noisy, W); 

figure(2);
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'interp','edgecolor', 'none');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
axis off;

%% Time-stepping for harmonic mapping.

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

    [u1, u2, u3] = cpSphere(unew1, unew2, unew3);
    
    t = kt*dt;
    error(kt) = norm([u1,u2,u3] - [u1_init, u2_init, u3_init]); 
end

%% plot difused mapping visualizing with the texture mapped image.
U_plot = [Eplot*u1, Eplot*u2, Eplot*u3];
[V, F, C] = tri_color(U_plot, W); 

figure(3);
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'interp','edgecolor', 'none');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
axis off;

%% plot the difference between the coordinates and original coordinates.
figure(4);
plot(1:numtimesteps,error,'k');
xlabel('$t$','Interpreter','latex','FontSize',20);
ylabel('$\|{\bf u}({\bf x},t)-{\bf u}_0({\bf x})\|_2$','Interpreter','latex','FontSize',20);

t_explicit = toc
