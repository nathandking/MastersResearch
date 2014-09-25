function [n,Qx,Qy,Qz,cpx,cpy,cpz,newx,newy]=CompGridCPThenMDS(dx)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% save results
save('dx0p1/Qx','Qx');
save('dx0p1/Qy','Qy');
save('dx0p1/Qz','Qz');
save('dx0p1/newx','newx');
save('dx0p1/newy','newy');
