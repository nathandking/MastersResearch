load('Sphere_MDSp05.mat');

%% use above flattened coordinates to place texture image on surface.
T = imread('../../Images/StJohns.jpg');
T = flipud(T);
T = fliplr(T);
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
% uv = [uu(:),vv(:)];
% newxy = [newx, newy];
% k = dsearchn(uv,newxy);
Ur = griddata(uu(:), vv(:), Tr, newxy(:,1), newxy(:,2),'natural');
Ug = griddata(uu(:), vv(:), Tg, newxy(:,1), newxy(:,2),'natural');
Ub = griddata(uu(:), vv(:), Tb, newxy(:,1), newxy(:,2),'natural');

Ur = Ur/max(Ur);
Ug = Ug/max(Ug);
Ub = Ub/max(Ub);

Urgb = zeros(size(cpX,1),3);
Urgb(sub_idx,:) = [Ur, Ug, Ub];

% figure(1);
% scatter3(cpX(:,1),cpX(:,2),cpX(:,3),20,Urgb,'fill');
% view([90 0]);
% axis([-1 1 -1 1 -1 1]);
% %colormap('gray');
% axis off

%% triangulation of surface using Delaunay triangulation.
[uni_cpX, ia, ~] = unique(cpX,'rows'); 
uni_Urgb = Urgb(ia,:);

DT = delaunayTriangulation(uni_cpX(:,1),uni_cpX(:,2),uni_cpX(:,3));
[Tri, Xb] = freeBoundary(DT);

[~,ib,ic]=intersect(uni_cpX, Xb, 'rows');
Tri_Urgb = uni_Urgb(ib,:);
TriP = Xb(ic,:);
iic(ic) = 1:length(ic);
Trin = iic(Tri);

% figure;
% trisurf(Trin,TriP(:,1),TriP(:,2),TriP(:,3),Tri_Urgb,'EdgeColor','none');
% view([90 0]);
% axis([-1 1 -1 1 -1 1]);
% axis off;

figure(1);
patch('Vertices', TriP, 'Faces', Trin, ... 
    'FaceVertexCData', Tri_Urgb, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');
view([90 0]);
axis([-1 1 -1 1 -1 1]);
axis off;