%% Find the grid point to center the interpolation scheme around.
function [Li,Lj]=FindIndex(dx,LL,cpx,cpy,Qx,Qy)

TRx=ceil((cpx+0.5*LL)/dx);
TRy=ceil((cpy+0.5*LL)/dx);

if cpx>=(Qx(TRx,TRy)-0.5*dx)
    Li=TRx;
else
    Li=TRx-1;
end

if cpy>=(Qy(TRx,TRy)-0.5*dx)
    Lj=TRy;
else
    Lj=TRy-1;
end

