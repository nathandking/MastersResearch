%% Definition of the closest point function for a circle.
function [cpx,cpy]=CP(x,y)

cpx=x./sqrt(x.^2+y.^2);
cpy=y./sqrt(x.^2+y.^2);

cpx(isnan(cpx))=1;
cpy(isnan(cpy))=0;