%% Level Set Method for Heat equation on a Sphere
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

% %% Banding: do calculation in a narrow band around the sphere
% dim = 3;    % dimension
% p = 3;      % interpolation degree
% order = 2;  % Laplacian order
% % "band" is a vector of the indices of the points in the computation
% % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% % the 1.0001 is a safety factor.
% bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
% band = find(abs(dist) <= bw*dx);
% 
% % store closest points in the band;
% cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
% xg = xx(band); yg = yy(band); zg = zz(band);


%% Function u in the embedding space
% this makes u into a vector, containing only points in the band
u1 = cpx;
u2 = cpy;
u3 = cpz;

%% Compute initial mapping of image onto sphere.
disp('Constructing initial map');
[newx, newy, xS, yS, zS, U] = InitialMap(41);
W = double(U);
disp('Initial map made');
%% Interpolation to get color onto computational points.
xS1 = xS(:); yS1 = yS(:); zS1 = zS(:);
size(u1)
Uc = griddata(xS1, yS1, zS1, W, u1(:), u2(:), u3(:),'nearest');
size(u1)
%% Visualize initial map.
% DT = delaunayTriangulation(u1,u2,u3);
% Tri = freeBoundary(DT);
% figure;
% trisurf(Tri,u1,u2,u3,Uc);

figure;
scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
colormap('copper');
axis([-1 1 -1 1 -1 1]);

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
figure;
scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
colormap('copper');
axis([-1 1 -1 1 -1 1]);

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
   

Tf = 0.05;
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
figure;
scatter3(u1(:),u2(:),u3(:),20,Uc,'fill');
colormap('copper');
axis([-1 1 -1 1 -1 1]);

toc