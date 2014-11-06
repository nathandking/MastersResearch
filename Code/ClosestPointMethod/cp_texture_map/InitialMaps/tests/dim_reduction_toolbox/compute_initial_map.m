clear all
close all

tic
n = 20;
method = 'Isomap';
[newx, newy, xS, yS, zS, U] = InitialMap(n, method);
W = double(U);
tfin = toc

save(strcat(method,'_StJohnsSphere_drtoolbox_',num2str(n+1),'.mat'),'newx','newy',...
    'xS','yS','zS','U','tfin');