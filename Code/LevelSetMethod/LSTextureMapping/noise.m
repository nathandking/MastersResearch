load('Uitest2.mat'); load('xtest2.mat'); load('ytest2.mat'); load('ztest2.mat');

Unoise = imnoise(Ui,'gaussian',1e-10);

figure;
s = surf(x,y,z,Unoise,'EdgeColor','none');