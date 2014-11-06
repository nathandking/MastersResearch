function [newxy] = MDS(M)
NumPts = size(M,1);
%% classical scaling algorithm
    J = eye(NumPts) - ones(NumPts)./(NumPts);
    B = -0.5 * J * M * J;
    % Find largest eigenvalues + their eigenvectors
    opts.tol = 1e-8;
    [Q, L,~] = eigs(B,2,'LM',opts);
    % Extract the coordinates.
    newx = sqrt(L(1,1)).*Q(:,1);
    newy = sqrt(L(2,2)).*Q(:,2);
    
    newxy = [newx, newy];