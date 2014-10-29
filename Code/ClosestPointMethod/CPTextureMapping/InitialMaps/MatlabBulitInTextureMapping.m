clear all;
close all;

StJohns = imread('StJohns.jpg');
Xhalf = [ones(1080,270,3),[ones(270,540,3); StJohns; ones(270,540,3)],ones(1080,270,3)];
[x,y,z] = sphere;
h = surf(x,y,z);
set(h, 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(Xhalf))
view([90 0])