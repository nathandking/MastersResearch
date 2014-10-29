function [not_c]= non_collinear(p)

not_c = 'true';
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
p4 = p(4,:);

if rank ([p1-p2 ; p3-p1]) < 2
    not_c = 'false';
end

if rank ([p4-p2 ; p3-p4]) < 2
    not_c = 'false';
end

if rank ([p1-p2 ; p4-p1]) < 2
    not_c = 'false';
end

if rank ([p1-p3 ; p4-p1]) < 2
    not_c = 'false';
end
