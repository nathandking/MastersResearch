function [patch, st_patch, cpX_tan_2d] = non_collinear_enclosing_patches(XS, cpX, st)

% compute all XS in tangent plane to cpX.
[XS_tan, cpX_tan] = Tp_sphere(XS,cpX);

XS_tan_2d = XS_tan(:,1:2);
cpX_tan_2d = cpX_tan(1:2)';

% check if 3 closest point are non collinear
if enclosing(XS_tan_2d(1:4,:), cpX_tan_2d)
    if non_collinear(XS_tan_2d(1:4,:))
        patch = XS_tan_2d(1:4,:);
        st_patch = st(1:4,:);
    else    % check the next closest point is non collinear instead of the fourth
       patch = [];
       st_patch = [];
       disp('not collinear')
    end
else     % check the next points until it is enclosing
    patch = [];
    st_patch = [];
end   