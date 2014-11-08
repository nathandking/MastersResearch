function save_initial_map()
tic
n = 200;
[newxy, sub_idx] = compute_initial_map(n);
tfin = toc;

save(strcat('Sphere_MDS_noncp',num2str(n+1),'.mat'),'newxy','sub_idx','tfin','n');%
quit;
end