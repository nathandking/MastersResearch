%% Level Set Method for a Circle
clear all;
tic

restoredefaultpath;
addpath('cp_matrices_Circle');
    %% construct the computational grid.
    dx=0.1;
    dy=dx;
    
    x = (-2.0:dx:2.0)';
    y = x;

    nx = length(x);
    ny = length(y);
    n = nx - 1;
    
    [Qx,Qy]=ndgrid(x,y);
    
    %% extend data using closest point extension.
    [cpx,cpy]=cpCircle(Qx,Qy); % x & y values of CPs to the computational grid.
    [th, r] = cart2pol(cpx,cpy); % define initial data on computational domain.
    u0 = cos(th); 
    u=u0;
    
    %% Construct an interpolation matrix for plotting on circle

% plotting grid on circle, using theta as a parameterization
thetas = linspace(0,2*pi,100)';
r = ones(size(thetas));
% plotting grid in Cartesian coords
[xp,yp] = pol2cart(thetas,r);
xp = xp(:); yp = yp(:);
Eplot = interp2_matrix(x, y, xp, yp,3,[],true);

figure(1); clf;
figure(2); clf;
figure(3); clf;

    
    %% compute the gradient of the level set function Phi.
    gPhix=Qx./sqrt(Qx.^2+Qy.^2);
    gPhix(isinf(1./gPhix))=0;
    gPhix(isnan(gPhix))=0;

    gPhiy=Qy./sqrt(Qx.^2+Qy.^2);
    gPhiy(isinf(1./gPhiy))=0;
    gPhiy(isnan(gPhiy))=0;

    % using signed distance function so Normal=grad(Phi).
    Nx=gPhix;
    Ny=gPhiy;
    
    %% solve PDE using level set method.
    vx=zeros(n+1);
    vy=zeros(n+1);
 
Tf = 1;
dt = 0.2*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

for kt = 1:numtimesteps
    
        % compute the gradient of u.
        vx(1:n,1:n)=(u(2:n+1,1:n)-u(1:n,1:n))/dx;
        vy(1:n,1:n)=(u(1:n,2:n+1)-u(1:n,1:n))/dx;
        % implement Neumann BCs.
        for i=1:n
            vx(i,n+1)=vx(i,n);
            vy(i,n+1)=vy(i,n);
        end
        for j=1:n
            vx(n+1,j)=vx(n,j);
            vy(n+1,j)=vy(n,j);   
        end
        vx(1,n+1)=vx(2,n); 
        vx(n+1,1)=vx(n,2);
        vx(n+1,n+1)=vx(n,n);
        vy(1,n+1)=vy(2,n); 
        vy(n+1,1)=vy(n,2);
        vy(n+1,n+1)=vy(n,n);
         
        % compute the projection.
        Projx=vx-(Nx.*vx+Ny.*vy).*Nx;
        Projy=vy-(Nx.*vx+Ny.*vy).*Ny;
        
        % compute the divergence of the projection.
        w(2:n+1,2:n+1)=(Projx(2:n+1,2:n+1)-Projx(1:n,2:n+1)...
            +Projy(2:n+1,2:n+1)-Projy(2:n+1,1:n))/dx;
        % implement Neumann BCs.
        for i=2:n+1
            w(i,1)=w(i,2);
        end
        for j=2:n+1
            w(1,j)=w(2,j);   
        end
        w(1,1)=w(2,2);
        w(1,n+1)=w(2,n); 
        w(n+1,1)=w(n,2);
    
        % do one step of forward Euler.
        u=u+dt*w;
        % implement Neumann BCs.
        for i=2:n
            u(i,1)=u(i,2);
            u(i,n+1)=u(i,n);
        end
        for j=2:n
            u(1,j)=u(2,j);
            u(n+1,j)=u(n,j);
        end
        u(1,1)=u(2,2);
        u(1,n+1)=u(2,n); 
        u(n+1,1)=u(n,2);
        u(n+1,n+1)=u(n,n);
        
        t = kt*dt;
  %%      
  if ( (kt < 10) || (mod(kt,10) == 0) || (kt == numtimesteps) )
    % plot over computation band
    plot2d_compdomain(u(:), cpx(:), cpy(:), dx, dx, 1)
    title( ['embedded domain: soln at time ' num2str(t) ...
            ', timestep #' num2str(kt)] );
    xlabel('x'); ylabel('y');
    %hold on
    plot(xp,yp,'k-', 'linewidth', 2);
    %axis equal;  axis tight;

    % plot value on circle
    set(0, 'CurrentFigure', 2);
    clf;
    circplot = Eplot*u(:);
    exactplot = exp(-t)*cos(thetas);
    plot(thetas, circplot);
    title( ['soln at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('u');
    hold on;
    % plot analytic result
    plot(thetas, exactplot, 'r--');
    legend('explicit Euler', 'exact answer', 'Location', 'SouthEast');
    error_circ_inf = max(abs( exactplot - circplot ))

    set(0, 'CurrentFigure', 3);
    clf;
    plot(thetas, circplot - exactplot);
    title( ['error at time ' num2str(t) ', on circle'] );
    xlabel('theta'); ylabel('error');

    %pause(0);
    drawnow();
  end
end
    
toc
