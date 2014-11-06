function [grad_signed] = grad_signed_Ellipsoid(x, y, z, AB)

grad = 2*[x/AB(1)^2, y/AB(2)^2, z];

grad_signed = zeros(size(x,1),3);
for i = 1:size(x,1)
    grad_signed(i,:) = grad(i,:)./norm(grad(i,:),2);
end

