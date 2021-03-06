tic
load('../flattening/Sphere_MDS_noncp201.mat');
%n = sqrt(size(XS,1))-1;

%% use above flattened coordinates to place texture image on surface.
T = imread('../../../Images/StJohns.jpg');
T = imresize(T, 0.25);
Tr = double(T(:,:,1));
Tg = double(T(:,:,2));
Tb = double(T(:,:,3));

Tr = Tr(:);
Tg = Tg(:);
Tb = Tb(:);

newx = newxy(:,1); newy = newxy(:,2);

% scale texture image to lie within (newx,newy) coordinates.
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

[uu,vv] = ndgrid(u,v);

%find indices of the closest points in (u,v) to (newx,newy).
Ur = griddata(uu(:), vv(:), Tr, newxy(:,1), newxy(:,2),'natural');
Ug = griddata(uu(:), vv(:), Tg, newxy(:,1), newxy(:,2),'natural');
Ub = griddata(uu(:), vv(:), Tb, newxy(:,1), newxy(:,2),'natural');

Ur = Ur/max(Ur);
Ug = Ug/max(Ug);
Ub = Ub/max(Ub);

Urgb = zeros((n+1)^2,3);
Urgb(sub_idx,:) = [Ur, Ug, Ub];
tfin_color = toc;
save(strcat('StJohns_Sphere_MDS_noncp_color',num2str(n+1),'.mat'),'Urgb','tfin_color');




