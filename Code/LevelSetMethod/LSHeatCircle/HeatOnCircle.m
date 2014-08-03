%% Level set method for heat equation on a circle
% This code implements the method in "Variational problems and partial
% differential equations on implicit surfaces", Bertalmio, Cheng, Osher and
% Sapiro. This code is just for the simple case when the surface S is a
% circle. To extend the original surface data into the embedding space, a
% constant normal extension is used.
clear all;
tic

%% loop over different dx values.
ddx=[0.2,0.1,0.05,0.025];

for q=1:2
    %% construct the computational grid.
    dx=ddx(q);
    dy=dx;
    dt=0.1*dx^2;
    LL=4;
    n=(LL/dx);
    x=linspace(-0.5*LL,0.5*LL,n+1);  y=linspace(-0.5*LL,0.5*LL,n+1);
    [Qx,Qy]=ndgrid(x,y);
    
    %% extend data using closest point extension.
    [cpx,cpy]=CP(Qx,Qy); % x & y values of CPs to the computational grid.
    u0=cpy./sqrt(cpx.^2+cpy.^2); % define initial data on computational domain.
    u=u0;
%     figure(1)
%     surf(Qx,Qy,u0)
    
    %% compute exact gradient of the level set function Phi.
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
    t=0:dt:1;
    for m=1:length(t)
    
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
        
        % plot surface at each time step.
%         figure(2)
%         surf(Qx,Qy,u)
%         hold off
%         pause(0.05)
    end
    
    %% interpolate solution to compare with exact solution.
    Li=zeros(n+1);
    Lj=zeros(n+1);

    for l=1:n
        for k=1:n
            % for each point this gives the grid for Lagrange interpolation.
            [Li(k,l),Lj(k,l)]=FindIndex(dx,LL,cpx(k,l),cpy(k,l),Qx,Qy);
        end
    end

    for l=1:n+1
        for k=1:n+1
            if Li(k,l)~=0
                Litmp=Li(k,l);
                Ljtmp=Lj(k,l);
            else
                continue;
            end
            
            Inx=find(Li-Litmp==0);
            Iny=find(Lj-Ljtmp==0);
            Ind=intersect(Inx,Iny);
            [row,col]=ind2sub([n+1,n+1],Ind);
            Xq=zeros(1,length(row));
            Yq=zeros(1,length(row));
        
            for a=1:length(row)
                Xq(a)=cpx(row(a),col(a));
                Yq(a)=cpy(row(a),col(a));
            end
            
            % points for Lagrange Interpolation.
            Lx=Qx(Litmp-2:Litmp+2,Ljtmp); 
            Ly=Qy(Litmp,Ljtmp-2:Ljtmp+2);  
            % solution values at the grid points.
            v=u(Litmp-2:Litmp+2,Ljtmp-2:Ljtmp+2);    
           
            uStmp=barylag2d(v,Lx,Ly,Xq,Yq);
            for a=1:length(row)
                u(row(a),col(a))=uStmp(a,a);
            end
            Li(Ind)=Li(Ind)-Litmp;
            Lj(Ind)=Lj(Ind)-Ljtmp;
        end
    end
    
    %% calculate error in numerical solution.
    error=zeros(n);
    for i=2:n
        for j=2:n
            [theta,rho]=cart2pol(x(i),y(j));
            uref(i,j)=exp(-t(end))*sin(theta);
            error(i,j)=abs(uref(i,j)-u(i,j));
        end
    end
    error=error(:);
    Error(q)=norm(error,inf)
end
% plot final numerical solution and error convergence.
figure(3)
surf(Qx(2:n,2:n),Qy(2:n,2:n),u(2:n,2:n))
hold on
surf(Qx(2:n,2:n),Qy(2:n,2:n),uref(2:n,2:n))

f1=figure(4);
loglog(ddx(1:q),Error)
xlabel('log(dx)','Interpreter','latex','FontSize',16);
ylabel('log(Error)','Interpreter','latex','FontSize',16);
axis([0.02 0.3 0.005 0.06]);
print(f1,'-deps','~/plot.eps');
hold on
plot([0.02 0.3],[0.005 0.06],'k--')
h=text(0.07,0.012,'Slope=1','FontSize',16,'Interpreter','latex');
set(h,'rotation',35)
order=log2(Error(1:end-1)./Error(2:end))'

toc
