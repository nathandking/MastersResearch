clear all
close all

M = 512;
N = 512;

u0=ones(M,N);

figure(1)
imagesc(u0); axis image; axis off; colormap(gray);
h1 = impoly;
h2 = imellipse;
h3 = imellipse;
h4 = impoly;
h5 = imellipse;
h6 = impoly;

umat = u0 - createMask(h1) - createMask(h2) + createMask(h3) -...
    createMask(h4) + createMask(h5) - createMask(h6);


figure(2)
imagesc(umat); axis image; axis off; colormap(gray);
umat(umat==-1) = 0;

imwrite(umat,'JackOLantern.jpg')