tic
load('Sphere_MDS_noncp301.mat');

%% use above flattened coordinates to place texture image on surface.
T = imread('../../Images/cboard.jpg');
T = double(T(:));

newx = newxy(:,1); newy = newxy(:,2);

% scale texture image to lie within (newx,newy) coordinates.
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

[uu,vv] = ndgrid(u,v);

%find indices of the closest points in (u,v) to (newx,newy).
X = [uu(:),vv(:)];
Y = [newxy(:,1), newxy(:,2)];
IDX = knnsearch(X,Y);
U = T(IDX);

%U = griddata(uu(:), vv(:), T, newxy(:,1), newxy(:,2),'nearest');
%U = U/max(U);

Ugray = zeros(size(XS,1),3);
Ugray(sub_idx) = U;
tfin_gray = toc;
save(strcat('cboard_Sphere_MDS_noncp_grayscale',num2str(sqrt(size(XS,1))),'.mat'),'Ugray','tfin_gray');