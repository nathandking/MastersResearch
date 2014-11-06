restoredefaultpath;
addpath('../cp_matrices_Sphere');

tic
dx = 0.05;
[newxy, cpX] = compute_initial_map(dx);
tfin = toc;

strdx = strrep(num2str(dx),'0.','p');
save(strcat('Sphere_MDS',strdx,'.mat'),'newxy','cpX','tfin');