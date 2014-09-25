%% First attempt to do manifold mapping with closest point method.
% This code does not implement the cp(x) on the plane first before the
% cp(x) onto the sphere. However since we are just doing the laplacian of
% this and the computational grid is uniform this should be the same result
% since the laplacian should not change points away from the boundary.
% However the general approach needs to be implemented.

restoredefaultpath;
addpath('cp_matrices_Sphere');
%%
% 3D example on a sphere
% Construct a grid in the embedding space

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
[xx yy zz] = meshgrid(x1d, y1d, z1d);
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
% u is a function defined on the grid (eg heat if solving the heat
% equation)

% assign some initial value (using initial value of cos (8*theta))
%[th, phi, r] = cart2sph(xx,yy,zz);
%u = cos(phi + pi/2);

% this makes u into a vector, containing only points in the band
u1 = cpxg;
u2 = cpyg;
u3 = cpzg;

%% Construct an interpolation matrix for closest point

% This creates a matrix which interpolates data from the grid x1d y1d,
% onto the points cpx cpy.

disp('Constructing interpolation and laplacian matrices');

E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);

% e.g., closest point extension:
%u = E*u;

%% Create Laplacian matrix for heat equation


L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);


%% Construct an interpolation matrix for plotting on sphere

% plotting grid on sphere, based on parameterization
% [xp,yp,zp] = sphere(64);
% xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
% [th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% % Eplot is a matrix which interpolations data onto the plotting grid
% Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);
% 
% figure(3); set(gcf,'Position', [410 700 800 800]);


%% Time-stepping for the heat equation

Tf = 0.02;
dt = 0.1*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
tic
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew1 = u1 + dt*(L*u1);
    unew2 = u2 + dt*(L*u2);
    unew3 = u3 + dt*(L*u3);
    
    % closest point extension
    u1 = E*unew1;
    u2 = E*unew2;
    u3 = E*unew3;
    t = kt*dt;
    
end

load('dx0p2/newx.mat');
load('dx0p2/newy.mat');
uvColor = checkerboard(10,2,2);
xmi = min(min(newx)); xma = max(max(newx));
ymi = min(min(newy)); yma = max(max(newy));

u = xmi:(xma-xmi)/(size(uvColor,1)-1):xma;
v = ymi:(yma-ymi)/(size(uvColor,2)-1):yma;

figure;
imagesc(uvColor);

%% relate (newx,newy) with corresponding colors on checkered image

Indx = round(0.5*(newx+abs(xmi))*(size(uvColor,1)-1)/xma)+1;
Indy = round(0.5*(newy+abs(ymi))*(size(uvColor,2)-1)/yma)+1;

for i = 1:length(Indx)
Ui(i) = uvColor(Indx(i),Indy(i));
end

DT = delaunayTriangulation(u1(:),u2(:),u3(:));
[F,V] = freeBoundary(DT);
figure;
patch('Faces',F,'Vertices',V);
%scatter3(u1(:),u2(:),u3(:),20,Ui(band))
t_explicit = toc
