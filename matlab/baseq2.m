function [fi,fix,fiy,fixx,fiyy,fixy,da]=baseq2(x,y,xc,yc)
%
%     As baseq, but also returning the second derivatives of the
%     quadratic basis functions.  They are constant on a straight-sided
%     triangle and are needed for the elementwise operator
%     div(d^3 grad v) in the stabilized Reynolds formulation.
%
%     xc,yc            element geometry node coordinates (3) for 1,3,6
%     da               area measure
%
%     6
%     | \
%     4   5
%     |    \
%     1--2--3
%

    da=polyarea(xc,yc);
    x1=xc(1);x2=xc(2);x3=xc(3);
    y1=yc(1);y2=yc(2);y3=yc(3);

    M=[1, x1, y1, x1*y1, x1^2, y1^2;
        1, (x1 + x2)/2,(y1 + y2)/2, ((x1 + x2)*(y1 + y2))/4, (x1 + x2)^2/4, (y1 + y2)^2/4;
        1, x2, y2, x2*y2, x2^2, y2^2;
        1, (x1 + x3)/2, (y1 + y3)/2, ((x1 + x3)*(y1 + y3))/4,(x1 + x3)^2/4, (y1 + y3)^2/4;
        1, (x2 + x3)/2,(y2 + y3)/2, ((x2 + x3)*(y2 + y3))/4, (x2 + x3)^2/4, (y2 + y3)^2/4;
        1, x3, y3, x3*y3, x3^2, y3^2];
    Mi=inv(M');

fi  =Mi*[1, x, y, x*y, x^2, y^2]';
fix =Mi*[0, 1, 0, y, 2*x, 0]';
fiy =Mi*[0, 0, 1, x, 0, 2*y]';
fixx=Mi*[0, 0, 0, 0, 2, 0]';
fiyy=Mi*[0, 0, 0, 0, 0, 2]';
fixy=Mi*[0, 0, 0, 1, 0, 0]';
