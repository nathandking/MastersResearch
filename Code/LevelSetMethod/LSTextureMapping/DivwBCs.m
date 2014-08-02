%% Computation of divergence of the projection into z = 0.
%--------------------------------------------------------------------------
%   * First order backward finite differences are used.
%   * Mirror boundary conditions are used here.
%   * Works only for arrays with the same size in all three dimensions.
%   * Works only for uniform stencils of the same size in three dimensions.
%--------------------------------------------------------------------------
function w = DivwBCs(Projx, Projy, Projz, dx)

% determine size of array.
n = size(Projx,1) - 1;
if (size(Projx,2) ~= n+1) || (size(Projx,3) ~= n+1)
    disp('problem with size of array');
end
%if (size(Projx) ~= size(Projy)) || (size(Projx) ~= size(Projz)) || (size(Projz) ~= size(Projy))
%    disp('problem with size of projection arrays');
%end

% backward difference for the divergence of the projection.
w(2:n+1,2:n+1,2:n+1)=(Projx(2:n+1,2:n+1,2:n+1)-Projx(1:n,2:n+1,2:n+1)...
            +Projy(2:n+1,2:n+1,2:n+1)-Projy(2:n+1,1:n,2:n+1)...
            +Projz(2:n+1,2:n+1,2:n+1)-Projz(2:n+1,2:n+1,1:n))/dx;
        
% implement Neumann BCs.
for i=2:n+1
    w(i,1,1) = w(i,2,2);
end
for j=2:n+1
    w(1,j,1) = w(2,j,2);   
end
for k=2:n+1
    w(1,1,k) = w(2,2,k);   
end
w(1,1,1) = w(2,2,2);
w(1,n+1,1) = w(2,n,2); 
w(n+1,1,1) = w(n,2,2);

