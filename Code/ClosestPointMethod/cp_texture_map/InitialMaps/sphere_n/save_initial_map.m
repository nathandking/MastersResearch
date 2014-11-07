tic
n = 400;
[newxy, XS, sub_idx] = compute_initial_map(n);
tfin = toc;

save(strcat('Sphere_MDS_noncp',num2str(n+1),'.mat'),'newxy','XS','sub_idx','tfin');