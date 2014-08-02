%% Computation of gradient of a three dimensional array.
%--------------------------------------------------------------------------
%   * First order forward finite differences are used.
%   * Mirror boundary conditions are used here.
%   * Works only for arrays with the same size in all three dimensions.
%   * Works only for uniform stencils of the same size in three dimensions.
%--------------------------------------------------------------------------
function [vx, vy, vz] = GradwBCs(u, dx)        

% determine size of array.
n = size(u,1) - 1;
if (size(u,2) ~= n+1) || (size(u,3) ~= n+1)
    disp('problem with size');
end

% forward difference for gradient of u.
vx(1:n,1:n,1:n) = (u(2:n+1,1:n,1:n) - u(1:n,1:n,1:n))/dx;
vy(1:n,1:n,1:n) = (u(1:n,2:n+1,1:n) - u(1:n,1:n,1:n))/dx;
vz(1:n,1:n,1:n) = (u(1:n,1:n,2:n+1) - u(1:n,1:n,1:n))/dx;
        
% implement Neumann BCs.
for i=1:n
    vx(i,n+1,n+1) = vx(i,n,n);
    vy(i,n+1,n+1) = vy(i,n,n);
    vz(i,n+1,n+1) = vz(i,n,n);
end
for j=1:n
    vx(n+1,j,n+1) = vx(n,j,n);
    vy(n+1,j,n+1) = vy(n,j,n);
    vz(n+1,j,n+1) = vz(n,j,n); 
end
for k = 1:n
    vx(n+1,n+1,k) = vx(n,n,k);
    vy(n+1,n+1,k) = vy(n,n,k);
    vz(n+1,n+1,k) = vz(n,n,k); 
end
vx(1,n+1,n+1) = vx(2,n,n);
vx(n+1,1,n+1) = vx(n,2,n);
vx(n+1,n+1,1) = vx(n,n,2);
vx(n+1,n+1,n+1) = vx(n,n,n);

vy(1,n+1,n+1) = vy(2,n,n); 
vy(n+1,1,n+1) = vy(n,2,n);
vy(n+1,n+1,1) = vy(n,n,2);
vy(n+1,n+1,n+1) = vy(n,n,n);
        
vz(1,n+1,n+1) = vz(2,n,n); 
vz(n+1,1,n+1) = vz(n,2,n);
vz(n+1,n+1,1) = vz(n,n,2);
vz(n+1,n+1,n+1) = vz(n,n,n);