% This is to test the MarchingCubes.m file to see if it can even replicate
% the surface of a sphere given its points.
dx = 0.2;
    %% construct the computational grid.
    dx = dx;
    dy = dx;
    dz = dx;
    
    x = (-2.0:dx:2.0)'; nx = length(x);
    y = x; 
    z = x; 
    n = nx-1;
    [Qx,Qy,Qz] = meshgrid(x,y,z);
    
    %% extend data using level set function for unit sphere.
    % Note: since a signed distance function is used, this is equivalent to
    % a closest point extension.
    [cpx,cpy,cpz,dist] = cpSphere(Qx,Qy,Qz);
    u1 = cpx;
    u2 = cpy;
    u3 = cpz;
    
    Ui = rand(nx,nx,nx);
% 
% [f,v,value] = MarchingCubes(Qx,Qy,Qz,sqrt(Qx.^2+Qy.^2+Qz.^2)-1,0,Ui);
% 
% 
% figure('color',[1 1 1])
% patch('Faces', f, 'Vertices', v, 'FaceVertexCData', value, ...
%           'FaceColor', 'flat', 'FaceLighting', 'phong', ...
%           'EdgeColor', 'none', 'userdata', value);

isosurface(Qx,Qy,Qz,sqrt(Qx.^2+Qy.^2+Qz.^2)-1,0,Ui);