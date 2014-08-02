%% Definition of the signed distance function for a circle.
function [SD]=SignDist(x,y,R)

[theta,rho]=cart2pol(x,y);

SD=rho-R;