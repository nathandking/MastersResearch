function [in] = enclosing(XS,cpX)

in = inpolygon(cpX(1),cpX(2),XS(:,1),XS(:,2));