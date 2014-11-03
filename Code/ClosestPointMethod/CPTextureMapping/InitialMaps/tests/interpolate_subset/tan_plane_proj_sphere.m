function [Tp, Y_3d] = tan_plane_proj_sphere(X,Y)

% compute the gradient of the signed distance function at Y.
gsd = grad_signed_dist_sphere(Y);

% project the points in X onto the tangent plane of Y.
Tp = zeros(size(X));
Y_3d = zeros(size(Y));
for i = 1:size(X,1)
Proj_Tp = eye(3) - gsd(i,:)'*gsd(i,:);
Tp(i,:) = Proj_Tp*X(i,:)';
Y_3d(i,:) = Proj_Tp*Y(i,:)';
end




