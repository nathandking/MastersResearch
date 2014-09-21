%% Level Set Method for Heat equation on a Sphere
tic

dx = 0.2;
[n,Qx,Qy,Qz,cpx,cpy,cpz,newx,newy] = CompGridCPThenMDS(dx);
%%
u1 = cpx;
u2 = cpy;
u3 = cpz;

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
    vx = zeros(n+1);
    vy = zeros(n+1);
    vz = zeros(n+1);
   

Tf = 0.1;
dt = 0.001*dx^2;
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create checkered image

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
Ui = reshape(Ui,21,21,21);

% figure;
% scatter3(cpx(:),cpy(:),cpz(:),20,Ui(:))
% 
% figure;
% scatter3(u1(:),u2(:),u3(:),20,Ui(:))
% figure;
% isosurface(Qx,Qy,Qz,sqrt(Qx.^2+Qy.^2+Qz.^2)-1,0,Ui);

figure;
isosurface(u1,u1,u3,sqrt(u1.^2+u2.^2+u3.^2)-1,0,Ui);

% [F,V,COLORS] = MarchingCubes(Qx,Qy,Qz,sqrt(Qx.^2+Qy.^2+Qz.^2)-1,0,Ui);


toc