%% compute geodesic distances of a sphere.
% takes ~7800 seconds with n = 120.
clear all;
close all;
tic
n = 40;
NumPts = (n+1)*(n+1);
% compute n+1 cartesian coordinates on a sphere.
[xTot, yTot, zTot] = sphere(4*(n-1));
x = xTot(1.5*n-1:2.5*n-1,1.5*n-1:2.5*n-1);
y = yTot(1.5*n-1:2.5*n-1,1.5*n-1:2.5*n-1);
z = zTot(1.5*n-1:2.5*n-1,1.5*n-1:2.5*n-1);

% convert cartesian coordinates to spherical coordinates.
%[a, e, r] = cart2sph(x,y,z);

% compute geodesic distances between n+1 points on sphere.
xx = x(:); yy = y(:); zz = z(:);
A = [xx, yy, zz];
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

Idx = reshape(Indx,n+1,n+1);
Idy = reshape(Indy,n+1,n+1);

for j = 1:n+1
    for i = 1:n+1
        Ui(i,j) = uvColor(Idx(i,j),Idy(i,j));
    end
end

figure;
imagesc(Ui)

%% put color of (newx,newy) on associated location on the sphere
figure;
s = surf(x,y,z,Ui,'EdgeColor','none');
% patch(surf2patch(s));
% delete(s)
% shading faceted; view(3)
toc