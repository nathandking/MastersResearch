%% Level Set Method for Heat equation on a Sphere
close all
clear all
tic
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

n = nx-1;
%% Find closest points on the surface
% For each point (x,y,z), we store the closest point on the sphere
% (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);
Qx = xx; Qy = yy; Qz = zz;


%% Function u in the embedding space
% this makes u into a vector, containing only points in the band

u1 = cpx;
u2 = cpy;
u3 = cpz;

%% load initial mapping of image onto sphere.

load('InitialMaps/SphereMDS101.mat');
W = double(U);

% Interpolation to get color onto computational points.
xS1 = xS(:); yS1 = yS(:); zS1 = zS(:);
Uc = griddata(xS1, yS1, zS1, W, u1(:), u2(:), u3(:),'nearest');

% Visualize initial map.
DT = delaunayTriangulation(u1(:),u2(:),u3(:));
Tri = freeBoundary(DT);
subplot(1,3,1);
trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
%axis off;

% subplot(1,3,1);
% scatter3(u3(:),u2(:),u1(:),20,Uc,'fill');
% axis([-1 1 -1 1 -1 1]);
% view([0 1 0]);
% colormap('copper');

%% Add noise to map.
N1 = 0.2*rand(length(u1(:)),1);
N2 = 0.2*rand(length(u2(:)),1);
N3 = 0.2*rand(length(u3(:)),1);
N1 = reshape(N1,size(u1));
N2 = reshape(N2,size(u2));
N3 = reshape(N3,size(u3));
u1 = u1 + N1;
u2 = u2 + N2;
u3 = u3 + N3;

% map noisy data back onto sphere.
[u1, u2, u3] = cpSphere(u1,u2,u3);

% visualize
DT = delaunayTriangulation(u1(:),u2(:),u3(:));
Tri = freeBoundary(DT);
subplot(1,3,2);
trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
%axis off;

% subplot(1,3,2);
% scatter3(u3(:),u2(:),u1(:),20,Uc,'fill');
% axis([-1 1 -1 1 -1 1]);
% view([0 1 0]);
% colormap('copper');

%% compute the gradient of the level set function Phi.
gPhix = Qx./sqrt(Qx.^2+Qy.^2+Qz.^2);
gPhix(isinf(1./gPhix)) = 0;
gPhix(isnan(gPhix)) = 0;

gPhiy = Qy./sqrt(Qx.^2+Qy.^2+Qz.^2);
gPhiy(isinf(1./gPhiy)) = 0;
gPhiy(isnan(gPhiy)) = 0;

gPhiz = Qz./sqrt(Qx.^2+Qy.^2+Qz.^2);
gPhiz(isinf(1./gPhiz)) = 0;
gPhiz(isnan(gPhiz)) = 0;

% using signed distance function so Normal=grad(Phi).
Nx = gPhix;
Ny = gPhiy;
Nz = gPhiy;
    
%% solve PDE using level set method.
Tf = 0.02;
dt = 0.001*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps
    %% compute gradient of u1, u2, u3.
    [v1x, v1y, v1z] = GradwBCs(u1,dx);
    [v2x, v2y, v2z] = GradwBCs(u2,dx);
    [v3x, v3y, v3z] = GradwBCs(u3,dx);
        
    %% divergence of the projections Proj1, Proj2, Proj3.
    w1 = DivwBCs(v1x, v1y, v1z, dx);
    w2 = DivwBCs(v2x, v2y, v2z, dx);
    w3 = DivwBCs(v3x, v3y, v3z, dx);
        
    %% project onto surface of sphere.
    uNormSq = u1.^2 + u2.^2 + u3.^2;
    SProj1 = w1 - (u1.^2.*w1 + u1.*u2.*w2 + u1.*u3.*w3)./uNormSq;
    SProj2 = w2 - (u2.*u1.*w1 + u2.^2.*w2 + u2.*u3.*w3)./uNormSq;
    SProj3 = w3 - (u3.*u1.*w1 + u3.*u2.*w2 + u3.^2.*w3)./uNormSq;
    
    %% one step of forward Euler.
    u1 = u1 + dt * SProj1;
    u2 = u2 + dt * SProj2;
    u3 = u3 + dt * SProj3;
        

        % implement Neumann BCs.
        % for u1.
        for i=2:n
            u1(i,1,1) = u1(i,2,2);
            u1(i,n+1,n+1) = u1(i,n,n);
        end
        for j=2:n
            u1(1,j,1) = u1(2,j,2);
            u1(n+1,j,n+1) = u1(n,j,n);
        end
        for k=2:n
            u1(1,1,k) = u1(2,2,k);
            u1(n+1,n+1,k) = u1(n,n,k);
        end
        u1(1,1,1) = u1(2,2,2);
        u1(1,n+1,n+1) = u1(2,n,n); 
        u1(n+1,1,n+1) = u1(n,2,n);
        u1(n+1,n+1,1) = u1(n,n,2);
        u1(n+1,n+1,n+1) = u1(n,n,n);
        % for u2.
        for i=2:n
            u2(i,1,1) = u2(i,2,2);
            u2(i,n+1,n+1) = u2(i,n,n);
        end
        for j=2:n
            u2(1,j,1) = u2(2,j,2);
            u2(n+1,j,n+1) = u2(n,j,n);
        end
        for k=2:n
            u2(1,1,k) = u2(2,2,k);
            u2(n+1,n+1,k) = u2(n,n,k);
        end
        u2(1,1,1) = u2(2,2,2);
        u2(1,n+1,n+1) = u2(2,n,n); 
        u2(n+1,1,n+1) = u2(n,2,n);
        u2(n+1,n+1,1) = u2(n,n,2);
        u2(n+1,n+1,n+1) = u2(n,n,n);
        % for u3.
        for i=2:n
            u3(i,1,1) = u3(i,2,2);
            u3(i,n+1,n+1) = u3(i,n,n);
        end
        for j=2:n
            u3(1,j,1) = u3(2,j,2);
            u3(n+1,j,n+1) = u3(n,j,n);
        end
        for k=2:n
            u3(1,1,k) = u3(2,2,k);
            u3(n+1,n+1,k) = u3(n,n,k);
        end
        u3(1,1,1) = u3(2,2,2);
        u3(1,n+1,n+1) = u3(2,n,n); 
        u3(n+1,1,n+1) = u3(n,2,n);
        u3(n+1,n+1,1) = u3(n,n,2);
        u3(n+1,n+1,n+1) = u3(n,n,n);
        
    t = kt*dt;
end

%% triangulation of surface using Delaunay triangulation.
DT = delaunayTriangulation(u1(:),u2(:),u3(:));
Tri = freeBoundary(DT);
subplot(1,3,3);
trisurf(Tri,u1(:),u2(:),u3(:),Uc,'EdgeColor','none');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
%axis off;

% subplot(1,3,3);
% scatter3(u3(:),u2(:),u1(:),20,Uc,'fill');
% view([0 1 0]);
% axis([-1 1 -1 1 -1 1]);
% colormap('copper');

toc