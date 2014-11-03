function [collinear]= test_not_collinearity(p1,p2,p3)

collinear = ones(size(p1,1),1);
for i = 1:size(p1,1)
    if rank ([p2(i,:)-p1(i,:) ; p3(i,:)-p1(i,:)]) < 2
        collinear(i) = 0;
    end
end

