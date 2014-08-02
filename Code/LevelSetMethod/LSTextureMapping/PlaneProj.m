%% Computation of projection into the plane z = 0.
%--------------------------------------------------------------------------
%   * Projects the gradient of three dimensional array into plane.
%--------------------------------------------------------------------------

function [Projx, Projy, Projz] = PlaneProj(vx, vy, vz)

   
Projx = vx;
Projy = vy;
Projz = vz - 1*vz*1;