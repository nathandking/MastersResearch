function [patch, st_patch, cpX_tan_2d] = non_collinear_enclosing_patches(XS, cpX, st)

% compute all XS in tangent plane to cpX.
[XS_tan, cpX_tan] = Tp_sphere(XS,cpX);

XS_tan_2d = XS_tan(:,1:2);
cpX_tan_2d = cpX_tan(1:2)';

%         figure(1)
%         axis equal
%         clf(1)
%         scatter(XS_tan_2d(:,1),XS_tan_2d(:,2),40,'fill')
%         hold on
%         scatter(XS_tan_2d(1:4,1),XS_tan_2d(1:4,2),40,'fill')
%         scatter(XS_tan_2d(5,1),XS_tan_2d(5,2),40,'fill')
%         %scatter(XS_tan_2d(6,1),XS_tan_2d(6,2),40,'fill')
%         scatter(cpX_tan_2d(1),cpX_tan_2d(2),40,'fill')
%        %pause(0.00005)



% check if 3 closest point are non collinear

if enclosing(XS_tan_2d(1:4,:), cpX_tan_2d) || enclosing(XS_tan_2d([1,3,2,4],:), cpX_tan_2d) || enclosing(XS_tan_2d([1,2,4,3],:), cpX_tan_2d)
    [not_c, id] = non_collinear(XS_tan_2d(1:4,:));
    if not_c
        patch = XS_tan_2d(1:4,:);
        st_patch = st(1:4,:);
    else    % check the next closest point is non collinear instead of the fourth
        patch = [];
        st_patch = [];
    end
else     % check the next points until it is enclosing
    
    % I first need to find 4 points that enclose.
    
    [not_c, id] = non_collinear(XS_tan_2d(1:4,:));
    if not_c
        next_try = 5;
        four_points = [XS_tan_2d(id,:); XS_tan_2d(next_try,:)];
      
        while (~enclosing(four_points, cpX_tan_2d) && ~enclosing(four_points([1,3,2,4],:), cpX_tan_2d) && ~enclosing(four_points([1,2,4,3],:), cpX_tan_2d)) || (~not_c)  
            next_try = next_try + 1;
            four_points = [XS_tan_2d(id,:); XS_tan_2d(next_try,:)]; 
            [not_c, ~] = non_collinear(four_points);
%         
%             figure(2)
%             clf(2)
%             scatter(four_points(:,1),four_points(:,2),40,'fill')
%             hold on
%             scatter(cpX_tan_2d(1),cpX_tan_2d(2),40,'fill')
%             pause(0.2)
            %enclosing(four_points, cpX_tan_2d)
%         
            if next_try >= size(XS_tan_2d,1)
                not_c;
                enclosing(four_points, cpX_tan_2d);
                error('not enough points in non enclosing case')
            end
        end
        patch = four_points;
        st_patch = [st(id,:); st(next_try,:)];
    else
        patch = [];
        st_patch = [];   
    end
end   