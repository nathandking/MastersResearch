clear all
close all

tic
n = 20;
[newx, newy, xS, yS, zS, U] = InitialMap(n);
W = double(U);
tfin = toc

save(strcat('StJohnsSphereMDS_subset',num2str(n+1),'.mat'),'newx','newy',...
    'xS','yS','zS','U','tfin');