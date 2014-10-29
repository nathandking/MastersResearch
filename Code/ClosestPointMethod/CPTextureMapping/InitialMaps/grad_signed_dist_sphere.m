function [gsd] = grad_signed_dist_sphere(Y)

% compute the gradient of the signed distance function at Y.
grad_signed_distx = Y(:,1)./sqrt(Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2);
grad_signed_disty = Y(:,2)./sqrt(Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2);
grad_signed_distz = Y(:,3)./sqrt(Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2);

gsd = [grad_signed_distx, grad_signed_disty, grad_signed_distz];