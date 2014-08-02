%% Level Set Method for Heat equation on a Sphere
clear all;
tic

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
    u1 = u0;
    u2 = u0;
    u3 = u0;
    
%% plotting grid on sphere, based on parameterization
[xp,yp,zp] = sphere(64);
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
[th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x, y, z, xp1, yp1, zp1);

figure(1); set(gcf,'Position', [410 700 800 800]);

    
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
   

Tf = 0.1;
dt = 0.0001*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps
tic
for kt = 1:numtimesteps
    
        %% compute gradient of u1, u2, u3.
        [v1x, v1y, v1z] = GradwBCs(u1,dx);
        [v2x, v2y, v2z] = GradwBCs(u2,dx);
        [v3x, v3y, v3z] = GradwBCs(u3,dx);
        
        %% project gradient of u1, u2, u3 into plane z = 0.
        [Proj1x, Proj1y, Proj1z] = PlaneProj(v1x,v1y,v1z);
        [Proj2x, Proj2y, Proj2z] = PlaneProj(v2x,v2y,v2z);
        [Proj3x, Proj3y, Proj3z] = PlaneProj(v3x,v3y,v3z);
        
        %% divergence of the projections Proj1, Proj2, Proj3.
        w1 = DivwBCs(Proj1x, Proj1y, Proj1z, dx);
        w2 = DivwBCs(Proj2x, Proj2y, Proj2z, dx);
        w3 = DivwBCs(Proj3x, Proj3y, Proj3z, dx);
        
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
        
    % plot value on sphere
    if (mod(kt,1000) == 0) || (kt < 10) || (kt == numtimesteps)
      figure(1);
      sphplot = Eplot*u1(:);
      
	  err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
      [t dt dx err]
      
      sphplot = reshape(sphplot, size(xp));
      surf(xp, yp, zp, sphplot);
      title( ['Level Set Method: soln at time ' num2str(t) ', kt= ' num2str(kt)] );
      xlabel('x'); ylabel('y'); zlabel('z');
      caxis([-1.05 1.05]);   % lock color axis
      axis equal; shading interp;
%      if ~exist OCTAVE_VERSION camlight left;
      colorbar;
      pause(0.001);
      end
   % end
        
    end
toc