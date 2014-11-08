function [newxy] = MDS(M)
%% classical scaling algorithm
    J = eye(size(M,1)) - ones(size(M,1))./(size(M,1));
    B = -0.5 * J * M * J;
    B = 0.5*(B+B');
    clear J;
    % Find largest eigenvalues + their eigenvectors
    opts.tol = 1e-8;
    [Q, L,~] = eigs(B,2,'LA',opts);
    clear B;
    % Extract the coordinates.
    newxy = [sqrt(L(1,1)).*Q(:,1), sqrt(L(2,2)).*Q(:,2)];
    clear L Q;
    
 