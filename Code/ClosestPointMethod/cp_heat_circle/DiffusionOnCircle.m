%% Closest Point Method for a Circle
% This code implements the method in "A simple embedding method for solving
% partial differential equations on surfaces", Ruuth and Merriman, 2007.
% This implementation is for the example in section 3.1 of the heat
% equation define on a circle.
clear all;
tic

% loop over different dx values.
ddx=[0.2,0.1,0.05,0.025,0.0125,0.00625];
for q=1:3
    %% Construct the computational grid.
    dx=ddx(q);
    dt=0.1*dx^2;
    LL=4;
    n=(LL/dx)-1;
    xx=linspace(-0.5*LL,0.5*LL,n+2);  yy=linspace(-0.5*LL,0.5*LL,n+2);
    x=xx(2:n+1);  y=yy(2:n+1);
    [Qx,Qy]=ndgrid(x,y);
    d=2;
    p=4;
    lambda=sqrt((d-1)*(0.5*(p+1))^2+(1+0.5*(p+1))^2)*dx;

    %% Extend initial surface data using the closest point function.
    [cpx,cpy]=CP(Qx,Qy); % x & y values of CPs to the computational grid.
    u=cpy./sqrt(cpx.^2+cpy.^2); % define initial data on computational domain.

    %%
    % Reorder into vector, with lexicographic (rowwise) ordering
    uvec=u(:);

    %%
    % Set up 2d discrete Laplacian in matrix form:
    onez = ones(n,1);
    D2 = spdiags(onez*[1 -2 1]/dx^2, -1:1, n, n);
    I = eye(n);
    A = kron(I,D2) + kron(D2,I);

    %%
    % Determine which points are in the computational band.
    Bd=sqrt((cpx-Qx).^2+(cpy-Qy).^2)-lambda<=0;

    %%
    t=0:dt:1;
    for m=1:length(t)
        % do one step of forward Euler, do not compute on the boundaries.
        uvec=uvec+dt*(A*uvec);
        u = reshape(uvec, n, n);
        % use fourth degree Lagrange interpolation to obtain approximate 
        % solution. Find the point closest to the point of interest and 
        % then choose two up and two down, two left and two right. This 
        % will make up the 5 by 5 grid.
        uS=u;
        Li=zeros(n);
        Lj=zeros(n);
        for l=1:n
            for k=1:n
                if Bd(k,l)==0
                    continue;
                else
                    % for each point on the circle this gives the grid for 
                    % Lagrange interpolation.
                    [Li(k,l),Lj(k,l)]=FindIndex(dx,LL,cpx(k,l),cpy(k,l),...
                        Qx,Qy);
                end
            end
        end
        
        for l=1:n
            for k=1:n
                if Bd(k,l)==0
                    continue;
                else
                    if Li(k,l)~=0
                        Litmp=Li(k,l);
                        Ljtmp=Lj(k,l);
                    else
                        continue;
                    end
                    
                    Inx=find(Li-Litmp==0);
                    Iny=find(Lj-Ljtmp==0);
                    Ind=intersect(Inx,Iny);
                    [row,col]=ind2sub([n,n],Ind);
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
                        uS(row(a),col(a))=uStmp(a,a);
                    end
                    
                    Li(Ind)=Li(Ind)-Litmp;
                    Lj(Ind)=Lj(Ind)-Ljtmp;
                end
            end
        end
        uvec=uS(:);
    end
    u = reshape(uvec, n, n);
    uref=exp(-t(end))*cpy./sqrt(cpx.^2+cpy.^2);
    error=zeros(n);
    for l=1:n
        for k=1:n
            if Bd(k,l)==0
                continue;
            else
                error(k,l)=abs(uref(k,l)-u(k,l));
            end
        end
    end
    error=error(:);
    Error(q)=norm(error,inf)
end
loglog(ddx(1:q),Error)
order=log2(Error(1:end-1)./Error(2:end))'
toc