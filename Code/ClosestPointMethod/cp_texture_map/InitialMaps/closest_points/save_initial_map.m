restoredefaultpath;
addpath('../../cp_matrices_Sphere');

tic
dx = 0.04;
[newxy, cpX, sub_idx] = compute_initial_map(dx);
tfin = toc;

strdx = strrep(num2str(dx),'0.','p');
save(strcat('Sphere_MDS',strdx,'.mat'),'newxy','cpX','sub_idx','tfin');