function [cpX_st] = projective_mapping2(cpX_tan_2d, p, st_p)

% implement projective mapping for quadrilateral-quadrilateral as outlines
% in Heckbert 1989.

M = zeros(8,8);
Msd = zeros(3,3);

M(1:4, 1:2) = [p(1,:); p(2,:); p(3,:); p(4,:)];
M(5:8, 4:5) = [p(1,:); p(2,:); p(3,:); p(4,:)];
M(1:4, 3) = ones(4,1);
M(5:8, 6) = ones(4,1);
M(:,7) = -[p(1,1)*st_p(1,1); p(2,1)*st_p(2,1); p(3,1)*st_p(3,1);...
    p(4,1)*st_p(4,1); p(1,1)*st_p(1,2); p(2,1)*st_p(2,2);...
    p(3,1)*st_p(3,2); p(4,1)*st_p(4,2)];
M(:,8) = -[p(1,2)*st_p(1,1); p(2,2)*st_p(2,1); p(3,2)*st_p(3,1);...
    p(4,2)*st_p(4,1); p(1,2)*st_p(1,2); p(2,2)*st_p(2,2);...
    p(3,2)*st_p(3,2); p(4,2)*st_p(4,2)];
   
b = [st_p(1,1); st_p(2,1); st_p(3,1); st_p(4,1); st_p(1,2); st_p(2,2); st_p(3,2); st_p(4,2)];
C = M\b;
   
Msd(:,1) = C(1:3);
Msd(:,2) = C(4:6);
Msd(:,3) = [C(7:8); 1];

cpX_st_temp = [cpX_tan_2d, 1]*Msd;
cpX_st = cpX_st_temp(1:2);
    


