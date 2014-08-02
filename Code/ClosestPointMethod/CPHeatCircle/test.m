dx=0.5; %dx=ddx(i);
dt=0.1*dx^2;
LL=3.6;
n=(LL/dx);
xx=linspace(-0.5*LL,0.5*LL,n+2);  yy=linspace(-0.5*LL,0.5*LL,n+2);
x=xx(2:n+1);  y=yy(2:n+1);
[Qx,Qy]=ndgrid(x,y)

%% Extend initial surface data using the closest point function.
[uSx,uSy]=CP(Qx,Qy) % x & y values of CPs to the computational grid.
u=uSy./sqrt(uSx.^2+uSy.^2); % define initial data on computational domain.