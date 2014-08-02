%% Closest Point Method for a Circle
% This code is a first implementation of a simple example for the Closest 
% Point Method.
clear all;

% Construct the computational grid.
%% Test construction of uniform mesh and matrix:
m = 4;
h = 1/(m+1);
x = linspace(0,1,m+2);              % Equivalent: x = 0:h:1;
y = linspace(0,1,m+2);
xi = x(2:m+1);                      % Interior points
yi = y(2:m+1);

%%
% Create mesh:
[xxi, yyi] = meshgrid(xi,yi);
disp('interior xx = '), disp(xxi)
disp('interior yy = '), disp(yyi)
figure(1), spy(xxi)                  % Display the grid

%%
% Reorder into vector, with rowwise ordering
%   (note: the colon operator : takes entries column by column):
xxv = xxi';     yyv = yyi';
xxv = xxv(:);   yyv = yyv(:);
disp('vectors are '),   disp([xxv yyv])

%%
% Set up 2d discrete Laplacian in matrix form:
onez = ones(m,1);
D2 = spdiags(onez*[1 -2 1], -1:1, m, m);
disp('1d differentiation matrix is '),      disp(full(D2))
I = eye(m);
disp('identity matrix is '),        disp(I)
L = kron(I,D2) + kron(D2,I);
disp('2d 5-point Laplacian matrix is '),    disp(full(L))
A = L/h^2;
figure(2), spy(A)

%% Now use larger m and solve Poisson's equation del^u = -1
m = 14;
h = 1/(m+1);
x = linspace(-3,3,m+2);              % Equivalent: x = 0:h:1;
y = linspace(-3,3,m+2);
xi = x(2:m+1);                      % Interior points
yi = y(2:m+1);

%%
% Create mesh, and reorder into vector:
[xxi, yyi] = meshgrid(xi,yi);
figure(3), spy(xxi)                 % Display the grid
xxv = xxi';     yyv = yyi';
xxv = xxv(:);   yyv = yyv(:);

%%
% Set up differentiation matrix:
onez = ones(m,1);
D2 = spdiags(onez*[1 -2 1], -1:1, m, m);
I = eye(m);
L = kron(I,D2) + kron(D2,I);
A = L/h^2;
figure(4), spy(A)

% Extend initial surface data using the closest point function.
[uSx,uSy]=CP(xxv,yyv); % x & y values of CPs to the computational grid.
u=sin(atan(uSy./uSx)); % define initial data on computational domain.

%%
% Solve linear system, and record the time:
tic
uvec = A\f;
t = toc;

%%
% Reshape vector of results onto 2d grid:
u = reshape(uvec, m, m);
uvec=u+dt*A*u;
uu = zeros(m+2, m+2);
uu(2:m+1, 2:m+1) = u1';              % Dirichlet BCs in rows/columns 1, m+2



[xx, yy] = meshgrid(x, y);

%%
% Plot solution:
figure(5), clf
mesh(xx, yy, uu)
title('Computed solution values')
xlabel('x')
ylabel('y')
zlabel('u')

figure(6), clf
clabel(contour(xx, yy, uu))
xlabel('x')
ylabel('y')
title(['Contours:  m = ' num2str(m) ', time for linear solver t = ' ...
    num2str(t)])
axis('square')
prism




%t=0:dt:1;
%for m=1:size(t)

% do one step of forward Euler, do not compute on the boundaries.
u(2:n-1,2:n-1)=u(2:n-1,2:n-1)+dt*(u(1:n-2,2:n-1)-2*u(2:n-1,2:n-1)+...
    u(3:n,2:n-1)+u(2:n-1,1:n-2)-2*u(2:n-1,2:n-1)+u(2:n-1,3:n))/dx^2;

% use fourth degree Lagrange interpolation to obtain approximate solution.

% Find the point closest to the point of interest and then choose two up
% and two down, two left and two right. This will make up the 5 by 5 grid.
uS=zeros(n);
for l=1:n
    for k=1:n
    % for each point on the circle this gives the grid for Lagrange interpolation.
        dtmp=IntSize;
        for j=1:n
            for i=1:n
                if sqrt((Qx(i,j)-uSx(k,l))^2+(Qy(i,j)-uSy(k,l))^2)<=dtmp
                dtmp=sqrt((Qx(i,j)-uSx(k,l))^2+(Qy(i,j)-uSy(k,l))^2);
                Li=i;
                Lj=j;
                end
            end
        end
        Lx=Qx(Li-2:Li+2,Lj-2:Lj+2);  % x points for Lagrange Interpolation.
        Ly=Qy(Li-2:Li+2,Lj-2:Lj+2);  % y points for Lagrange Interpolation.
        v=u(Li-2:Li+2,Lj-2:Lj+2);    % solution values at the grid points.

        % figure(1)
        % plot(Lx,Ly,'*')
        % hold on
        % plot(uSx(k,l),uSy(k,l),'d')
        
        % Create five interpolating polynomials along x with a fixed y. At
        % the x-coordinate of the point on the circle, compute the solution 
        % values using these five interpolating polynomials. Then create
        % one interpolating polynomial along y to obtain the solution value
        % at the point on the circle.
        yi=zeros(1,5);
        for i=1:5
        yi(i)=polyinterp(Lx(1,:),v(i,:),uSx(k,l));
        end
        uS(k,l)=polyinterp(Ly(:,1),yi,uSy(k,l));
    end
end
%m=m+1;
%end

uref=exp(-dt)*sin(atan(uSy./uSx));
Error=norm(abs(uS(2:n-1,2:n-1)-uref(2:n-1,2:n-1)),inf);
OutError=norm(abs(uS(2:n-1,[2:floor(0.5*n)-2 ceil(0.5*n)+2:n-1])-...
    uref(2:n-1,[2:floor(0.5*n)-2 ceil(0.5*n)+2:n-1])),inf);