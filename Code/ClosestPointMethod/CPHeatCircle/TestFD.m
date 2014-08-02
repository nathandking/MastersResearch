
clear all;
tic

% loop over different dx values.
ddx=[0.2,0.1,0.05,0.025,0.0125,0.00625];
for i=1:6

%% Construct the computational grid.
dx=ddx(i);
dt=0.1*dx^2;
LL=2;
n=(LL/dx)-1;
xx=linspace(0,LL,n+2);  yy=linspace(0,LL,n+2);
x=xx(2:n+1);  y=yy(2:n+1);
[Qx,Qy]=ndgrid(x,y);

%% Initial data.
u0=sin(Qx*pi/LL).*sin(Qy*pi/LL);
%surf(x,y,u0)
%%
% Reorder into vector, with lexicographic (rowwise) ordering
uvec=u0(:);

%%
% Set up 2d discrete Laplacian in matrix form:
onez = ones(n,1);
D2 = spdiags(onez*[1 -2 1], -1:1, n, n);
I = eye(n);
L = kron(I,D2) + kron(D2,I);
A = L/dx^2;

%%

t=0:dt:0.1;
for m=1:length(t)
% do one step of forward Euler, do not compute on the boundaries.
uvec=uvec+dt*(A*uvec);
end
u=reshape(uvec,n,n);

lambda2=2*(pi/LL)^2;
uref=exp(-lambda2*t(end))*sin(Qx*pi/LL).*sin(Qy*pi/LL);
%Error=norm((uu(2:n-1,2:n-1)-uref(2:n-1,2:n-1))./uref(2:n-1,2:n-1),inf);
%Error=max(max(abs((uS(2:n-1,2:n-1)-uref(2:n-1,2:n-1))./uref(2:n-1,2:n-1))));

Error(i)=max(max(abs((u(2:n-1,2:n-1)-uref(2:n-1,2:n-1))./uref(2:n-1,2:n-1))));
end
figure(1)
surf(x,y,uref)
figure(2)
loglog(ddx,Error)
order=log2(Error(1:5)./Error(2:6))'
toc

