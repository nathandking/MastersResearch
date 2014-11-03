function [Tp, cpX_tan] = Tp_sphere(XS,cpX)

% compute the gradient of the signed distance function at Y.
gsd = grad_signed_dist_sphere(cpX);

% project the points in X onto the tangent plane of Y.
Tp = zeros(size(XS));
cpX_tan = zeros(size(cpX));
for i = 1:size(XS,1)
    Proj_Tp = eye(3) - gsd'*gsd;
    Tp(i,:) = Proj_Tp*XS(i,:)';
    cpX_tan = Proj_Tp*cpX';
end