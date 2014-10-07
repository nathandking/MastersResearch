clear all
close all

tic
n = 100;
[newx, newy, xS, yS, zS, U] = InitialMap(n);
W = double(U);
tfin = toc

save(strcat('InitialMaps/SphereMDS',num2str(n+1),'.mat'),'newx','newy',...
    'xS','yS','zS','U','tfin');