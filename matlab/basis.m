function [fi,fix,fiy,area]=basis(x,y,xc,yc)
% Scalar evaluation of basis functions

x2=xc(2)-xc(1);
y2=yc(2)-yc(1);
x3=xc(3)-xc(1);
y3=yc(3)-yc(1);
detj=x2*y3-x3*y2;
xr=x-xc(1);
yr=y-yc(1);

fi=[0;xr*y3-x3*yr;x2*yr-xr*y2]/detj;
fi(1)=1-fi(2)-fi(3);

fix=[y2-y3;y3;-y2]/detj;

fiy=[x3-x2;-x3;x2]/detj;


area=(detj/2);
