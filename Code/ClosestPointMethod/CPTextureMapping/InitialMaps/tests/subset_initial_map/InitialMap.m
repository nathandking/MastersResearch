function [newx, newy, xS, yS, zS, U] = InitialMap(n)
%% Define coordinates on surface.
[xS,yS,zS] = sphere(n);

%% Compute geodesic distance between pairs of points.
% This below only works for a sphere!
S = [xS(:), yS(:), zS(:)];
[XS,ia,ic] = unique(S,'rows'); 
A = XS;
NumPts = size(XS,1);
M = zeros(NumPts);
for j = 1:NumPts
    for i = 1:NumPts
        M(i,j) = (real(acos(sum(A(i,:,:).*A(j,:,:))))).^2;
    end
end


    J = eye(NumPts) - ones(NumPts)./(NumPts);
    B = -0.5 * J * M * J;

%% Apply multidimensional scaling to determine flattened coordinates.
[newx, newy] = MDS(B);

%% use above flattened coordinates to place texture image on surface.
T = imread('../../../Images/StJohns.jpg');
T = rgb2gray(T);
Tg = double(T(:));


% scale texture image to lie within (newx,newy) coordinates.
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

[uu,vv] = ndgrid(u,v);

%find indices of the closest points in (u,v) to (newx,newy).
% uv = [uu(:),vv(:)];
% newxy = [newx, newy];
% k = dsearchn(uv,newxy);
U = griddata(uu(:), vv(:), Tg, newx, newy,'natural');
%U = Tg(k);
