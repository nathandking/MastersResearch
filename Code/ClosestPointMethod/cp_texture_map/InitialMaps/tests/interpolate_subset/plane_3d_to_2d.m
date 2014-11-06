function [p1_2d, p2_2d, p3_2d, p4_2d, Y_2d] = plane_3d_to_2d(p1,p2,p3,p4,Y)

% perform a change of basis to obtain a 2d plane from the 3d plane.
N = grad_signed_dist_sphere(Y);

W = p1 - p2;

% allocate memory.
U = zeros(size(Y,1),3);
V = zeros(size(Y,1),3);
p1_2d = zeros(size(Y,1),2);
p2_2d = zeros(size(Y,1),2);
p3_2d = zeros(size(Y,1),2);
p4_2d = zeros(size(Y,1),2);
Y_2d = zeros(size(Y,1),2);

for i = 1:size(Y,1)
U(i,:) = W(i,:) - (W(i,:)*N(i,:)')*N(i,:);
U(i,:) = U(i,:)./norm(U(i,:),2);
V(i,:) = cross(N(i,:)',U(i,:)');
p1_2d(i,:) = [p1(i,:)*U(i,:)', p1(i,:)*V(i,:)'];
p2_2d(i,:) = [p2(i,:)*U(i,:)', p2(i,:)*V(i,:)'];
p3_2d(i,:) = [p3(i,:)*U(i,:)', p3(i,:)*V(i,:)'];
p4_2d(i,:) = [p4(i,:)*U(i,:)', p4(i,:)*V(i,:)'];
Y_2d(i,:) = [Y(i,:)*U(i,:)', Y(i,:)*V(i,:)'];
end


