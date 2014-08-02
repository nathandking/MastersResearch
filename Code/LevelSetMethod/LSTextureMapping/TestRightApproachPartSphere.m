%% Level Set Method for Heat equation on a Sphere
clear all;
tic

    %% construct the computational grid.
    dx = 0.5;
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
 
    [azimuth,elevation,radius] = cart2sph(cpx,cpy,cpz);
    aG = find(azimuth(:) >= 0);  
    aL = find(azimuth(:) <= pi/2); 
    eG = find(elevation(:) >= -pi/4);
    eL = find(elevation(:) <= pi/4);
    S = [length(aG), length(aL), length(eG), length(eL)];
    [minS, minLoc] = min(S);
    
    if minLoc == 1
        minarray = aG;
        a1 = aL;
        a2 = eG;
        a3 = eL;
    elseif minLoc == 2
        minarray = aL;
        a1 = aG;
        a2 = eG;
        a3 = eL;
    elseif minLoc == 3
        minarray = eG;
        a1 = aG;
        a2 = aL;
        a3 = eL;
    else
        minarray = eL;
        a1 = aG;
        a2 = aL;
        a3 = eG;
    end
    
    [Lia1,Loca1] = ismember(minarray,a1);
    minarray1 = a1(Loca1(Loca1~=0));
    [Lia2,Loca2] = ismember(minarray1,a2);
    minarray2 = a2(Loca2(Loca2~=0));
    [Lia3,Loca3] = ismember(minarray2,a3);
    MinArray = a3(Loca3(Loca3~=0));
    
    Cx = cpx(MinArray);
    Cy = cpy(MinArray);
    Cz = cpz(MinArray);
    
    u1 = reshape(Cx,length(MinArray),9);
    u2 = reshape(Cy,length(MinArray),9);
    u3 = reshape(Cz,length(MinArray),9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %% compute geodesic distances of a sphere.

% convert cartesian coordinates to spherical coordinates.
%[a, e, r] = cart2sph(x,y,z);

% compute geodesic distances between n+1 points on sphere.
xx = u1(:); yy = u2(:); zz = u3(:);
A = [xx, yy, zz];
NumPts = length(xx);
M = zeros(NumPts);
for j = 1:NumPts
    for i = 1:NumPts
        M(i,j) = (real(acos(sum(A(i,:,:).*A(j,:,:))))).^2;
    end
end

%% classical scaling algorithm
    J = eye(NumPts) - ones(NumPts)./(NumPts);
    B = -0.5 * J * M * J;
    % Find largest eigenvalues + their eigenvectors
    [Q, L,flag] = eigs(B,2,'LM');
    % Extract the coordinates.
    newx = sqrt(L(1,1)).*Q(:,1);
    newy = sqrt(L(2,2)).*Q(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

uvColor = checkerboard;
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
% Idx = reshape(Indx,sqrt(NumPts),sqrt(NumPts));
% Idy = reshape(Indy,sqrt(NumPts),sqrt(NumPts));

% for j = 1:n+1
%     for i = 1:n+1
%         Ui(i,j) = uvColor(Idx(i,j),Idy(i,j));
%     end
% end

% figure;
% imagesc(Ui)
figure;
scatter3(u1(:),u2(:),u3(:),1,Ui(:))

figure;
fill3(u1(:),u2(:),u3(:),Ui(:))
%% put color of (newx,newy) on associated location on the sphere
% figure;
% surf(cpx,cpy,cpz,Ui,'EdgeColor','none');
% 
% figure;
% surf(u1,u2,u3,Ui,'EdgeColor','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc