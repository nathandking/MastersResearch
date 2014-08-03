%% Level set method for heat equation on a sphere.
% This code implements the method in "Variational problems and partial
% differential equations on implicit surfaces", Bertalmio, Cheng, Osher and
% Sapiro. This code is just for the simple case when the surface S is a
% sphere. To extend the original surface data into the embedding space, a
% constant normal extension is used. 
%
% This code uses the interpolation matrices from cp_matrices and plot 
% figures similar to example_heat_circle.m by Steve Ruuth and Colin 
% Macdonald.
clear all;
tic

restoredefaultpath;
addpath('cp_matrices_Sphere');
%% construct the computational grid.
dx = 0.2;
dy = dx;
dz = dx;
dt = 0.1*dx^2;
    
x = (-2.0:dx:2.0)'; nx = length(x);
y = x; ny = length(y);
z = x; nz = length(z);
n = nx-1;
[Qx,Qy,Qz] = ndgrid(x,y,z);
    
%% extend data using level set function for unit sphere.
% Note: since a signed distance function is used, this is equivalent to
% a closest point extension.
[cpx,cpy,cpz,dist] = cpSphere(Qx,Qy,Qz);
    
% define initial data on computational domain.
[theta,eta,r] = cart2sph(cpx,cpy,cpz);
u0 = cos(eta+pi/2); 
u = u0;
    
%% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x, y, z, xp1, yp1, zp1);

figure(2); set(gcf,'Position', [410 700 800 800]);
    
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
vx = zeros(nx);
vy = zeros(ny);
vz = zeros(nz);
    
Tf = 0.01;
dt = 0.0001*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
tic
for kt = 1:numtimesteps
    
    %% forward difference for gradient of u.
    vx(1:n,1:n,1:n) = (u(2:n+1,1:n,1:n) - u(1:n,1:n,1:n))/dx;
    vy(1:n,1:n,1:n) = (u(1:n,2:n+1,1:n) - u(1:n,1:n,1:n))/dy;
    vz(1:n,1:n,1:n) = (u(1:n,1:n,2:n+1) - u(1:n,1:n,1:n))/dz;
    % implement Neumann BCs.
    for i=1:n
        vx(i,n+1,n+1) = vx(i,n,n);
        vy(i,n+1,n+1) = vy(i,n,n);
        vz(i,n+1,n+1) = vz(i,n,n);
    end
    for j=1:n
        vx(n+1,j,n+1) = vx(n,j,n);
        vy(n+1,j,n+1) = vy(n,j,n);
        vz(n+1,j,n+1) = vz(n,j,n); 
    end
    for k = 1:n
        vx(n+1,n+1,k) = vx(n,n,k);
        vy(n+1,n+1,k) = vy(n,n,k);
        vz(n+1,n+1,k) = vz(n,n,k); 
    end
    vx(1,n+1,n+1) = vx(2,n,n);
    vx(n+1,1,n+1) = vx(n,2,n);
    vx(n+1,n+1,1) = vx(n,n,2);
    vx(n+1,n+1,n+1) = vx(n,n,n);
    
    vy(1,n+1,n+1) = vy(2,n,n); 
    vy(n+1,1,n+1) = vy(n,2,n);
    vy(n+1,n+1,1) = vy(n,n,2);
    vy(n+1,n+1,n+1) = vy(n,n,n);
        
    vz(1,n+1,n+1) = vz(2,n,n); 
    vz(n+1,1,n+1) = vz(n,2,n);
    vz(n+1,n+1,1) = vz(n,n,2);
    vz(n+1,n+1,n+1) = vz(n,n,n);
         
    %% project onto surface.
    Projx = vx - (Nx.*vx + Ny.*vy + Nz.*vz).*Nx;
    Projy = vy - (Nx.*vx + Ny.*vy + Nz.*vz).*Ny;
    Projz = vz - (Nx.*vx + Ny.*vy + Nz.*vz).*Nz;
        
    %% backward difference for the divergence of the projection.
    w(2:n+1,2:n+1,2:n+1)=(Projx(2:n+1,2:n+1,2:n+1)-Projx(1:n,2:n+1,2:n+1)...
        +Projy(2:n+1,2:n+1,2:n+1)-Projy(2:n+1,1:n,2:n+1)...
        +Projz(2:n+1,2:n+1,2:n+1)-Projz(2:n+1,2:n+1,1:n))/dx;
    % implement Neumann BCs.
    for i=2:n+1
        w(i,1,1) = w(i,2,2);
    end
    for j=2:n+1
        w(1,j,1) = w(2,j,2);   
    end
    for k=2:n+1
        w(1,1,k) = w(2,2,k);   
    end
    w(1,1,1) = w(2,2,2);
    w(1,n+1,1) = w(2,n,2); 
    w(n+1,1,1) = w(n,2,2);
    
    %% one step of forward Euler.
    u=u+dt*w;
    % implement Neumann BCs.
    for i=2:n
        u(i,1,1) = u(i,2,2);
        u(i,n+1,n+1) = u(i,n,n);
    end
    for j=2:n
        u(1,j,1) = u(2,j,2);
        u(n+1,j,n+1) = u(n,j,n);
    end
    for k=2:n
        u(1,1,k) = u(2,2,k);
        u(n+1,n+1,k) = u(n,n,k);
    end
    u(1,1,1) = u(2,2,2);
    u(1,n+1,n+1) = u(2,n,n); 
    u(n+1,1,n+1) = u(n,2,n);
    u(n+1,n+1,1) = u(n,n,2);
    u(n+1,n+1,n+1) = u(n,n,n);
        
    t = kt*dt;
        
    %% plot value on sphere
    if (mod(kt,100) == 0) || (kt < 10) || (kt == numtimesteps)
        figure(2);
        sphplot = Eplot*u(:);
     
        err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) /...
            norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
        [t dt dx err]
      
        sphplot = reshape(sphplot, size(xp));
        surf(xp, yp, zp, sphplot);
        title( ['Level Set Method: soln at time ' num2str(t) ', kt= ' num2str(kt)] );
        xlabel('x'); ylabel('y'); zlabel('z');
        caxis([-1.05 1.05]);   % lock color axis
        axis equal; shading interp;
        colorbar;
        pause(0.001);
    end
end
toc