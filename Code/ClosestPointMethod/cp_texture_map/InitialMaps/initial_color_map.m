load('Sphere_MDSp2.mat');

%% use above flattened coordinates to place texture image on surface.
T = imread('../Images/BlurredThreeShapes512.jpg');
%T = rgb2gray(T);
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
U = griddata(uu(:), vv(:), Tg, newxy(:,1), newxy(:,2),'natural');
%U = Tg(k);

figure;
scatter3(cpX(:,1),cpX(:,2),cpX(:,3),20,U,'fill');
view([0 1 0]);
axis([-1 1 -1 1 -1 1]);
colormap('gray');
axis off