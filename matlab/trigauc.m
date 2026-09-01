function [gcx,gcy,gv]=trigauc(xc,yc,poldegree)
% function [gc,gv]=trigauc(xc,yc,poldegree)
%
%     calculate coordinates and weight for Gauss quadrature
%      in the triangle with vertices in xc(*),yc(*)
%
%     poldegree: degree of polynomial to be integrated exactly (1,2,3,4,5,7,8 or 11)
%
switch(poldegree)
case 1
	gv=1;
	gcx=(xc(1)+xc(2)+xc(3))/3;
	gcy=(yc(1)+yc(2)+yc(3))/3;
case 2
	gv=[1/3,1/3,1/3];
	gcx=[(xc(1)+xc(2))/2,(xc(2)+xc(3))/2,(xc(3)+xc(1))/2];
	gcy=[(yc(1)+yc(2))/2,(yc(2)+yc(3))/2,(yc(3)+yc(1))/2];
case 3
	gv=[3,3,3,8,8,8,27]/60;
	gcx=[xc(1),xc(2),xc(3),(xc(1)+xc(2))/2,(xc(2)+xc(3))/2, ...
	(xc(3)+xc(1))/2,(xc(1)+xc(2)+xc(3))/3];
	gcy=[yc(1),yc(2),yc(3),(yc(1)+yc(2))/2,(yc(2)+yc(3))/2, ...
	(yc(3)+yc(1))/2,(yc(1)+yc(2)+yc(3))/3];
case 4
	gcx=zeros(1,6);gcy=gcx;gv=gcx;
	x1=0.8168475729804585;y1=0.09157621350977073;z1=1-x1-y1;  % multiplicity 3
	w1=0.3298552309659655/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(1:3)=xv;gcy(1:3)=yv;gv(1:3)=vv;
	x1=0.1081030181680702;y1=0.4459484909159649;z2=1-x1-y1;  % multiplicity 3
 	w1=0.6701447690340345/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(4:6)=xv;gcy(4:6)=yv;gv(4:6)=vv;
case 5
	x1=0.101286507323456;
	x2=0.797426985353087;
	x3=x1;
	x4=0.470142064105115;
	x5=x4;
	x6=0.059715871789770;
	x7=1/3;
	xsi=[x1,x2,x3,x4,x5,x6,x7];eta=[x1,x1,x2,x6,x4,x4,x7];
	w1=0.125939180544827;w4=0.132394152788506;w7=0.225;
	gv=[w1,w1,w1,w4,w4,w4,w7];
    gcx=[xsi*xc(1)+eta*xc(2)+(1-xsi-eta)*xc(3)];
    gcy=[xsi*yc(1)+eta*yc(2)+(1-xsi-eta)*yc(3)];
case 6
	gcx=zeros(1,12);gcy=gcx;gv=gcx;
	x1=0.5014265096581342;y1=0.2492867451709329;z1=1-x1-y1;  % multiplicity 3
	w1=0.3503588271790222/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(1:3)=xv;gcy(1:3)=yv;gv(1:3)=vv;
	x2=0.8738219710169965;y2=0.06308901449150177;z2=1-x2-y2;  % multiplicity 3
 	w2=0.1525347191106164/3;
	[xv,yv,vv]=multiplicity(x2,y2,w2,xc,yc,3);
	gcx(4:6)=xv;gcy(4:6)=yv;gv(4:6)=vv;
	x3=0.6365024991213939;y3=0.05314504984483216;z3=1-x3-y3;  % multiplicity 6
	w3=0.4971064537103575/6;
	[xv,yv,vv]=multiplicity(x3,y3,w3,xc,yc,6);
	gcx(7:12)=xv;gcy(7:12)=yv;gv(7:12)=vv;
case 7
	x1=0.065130102902216;
	x2=0.869739794195568;
	x3=x1;
	x4=0.312865496004874;
	x5=0.638444188569810;
	x6=0.048690315425316;
	x7=x5;
	x8=x4;
	x9=x6;
	x10=0.260345966079040;
	x11=0.479308067841920;
	x12=x10;
	x13=1/3;
	xsi=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13];
	eta=[x1,x1,x2,x6,x4,x5,x6,x5,x4,x10,x10,x11,x13];
	w1=0.053347235608838;w4=0.077113760890257;w10=0.175615257433208;
	w13=-0.149570044467682;
	gv=[w1,w1,w1,w4,w4,w4,w4,w4,w4,w10,w10,w10,w13];
    gcx=[xsi*xc(1)+eta*xc(2)+(1-xsi-eta)*xc(3)];
    gcy=[xsi*yc(1)+eta*yc(2)+(1-xsi-eta)*yc(3)];
case 8
	gcx=zeros(1,16);gcy=gcx;gv=gcx;
	x1=1/3;y1=1/3;z1=1-x1-y1;  % multiplicity 1
	w1=0.1443156076777862;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,1);
	gcx(1)=xv;gcy(1)=yv;gv(1)=vv;

	x1=0.08141482341455413;y1=0.4592925882927229;z1=1-x1-y1;  % multiplicity 3
	w1=0.2852749028018549/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(2:4)=xv;gcy(2:4)=yv;gv(2:4)=vv;

	x1=0.8989055433659379;y1=0.05054722831703103;z1=1-x1-y1;  % multiplicity 3
 	w1=0.09737549286959440/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(5:7)=xv;gcy(5:7)=yv;gv(5:7)=vv;

	x1=0.6588613844964797;y1=0.1705693077517601;z1=1-x1-y1;  % multiplicity 3
 	w1=0.3096521116041552/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(8:10)=xv;gcy(8:10)=yv;gv(8:10)=vv;

	x1=0.0083947774099572110;y1=0.7284923929554041;z1=1-x1-y1;  % multiplicity 6
	w1=0.1633818850466092/6;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,6);
	gcx(11:16)=xv;gcy(11:16)=yv;gv(11:16)=vv;
case 11
	gcx=zeros(1,28);gcy=gcx;gv=gcx;
	x1=0.9480217181434233;y1=0.02598914092828833;z1=1-x1-y1;  % multiplicity 3
	w1=0.02623293466120857/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(1:3)=xv;gcy(1:3)=yv;gv(1:3)=vv;

	x1=0.8114249947041546;y1=0.09428750264792270;z1=1-x1-y1;  % multiplicity 3
 	w1=0.1142447159818060/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(4:6)=xv;gcy(4:6)=yv;gv(4:6)=vv;

	x1=0.01072644996557060;y1=0.4946367750172147;z1=1-x1-y1;  % multiplicity 3
	w1=0.05656634416839376/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(7:9)=xv;gcy(7:9)=yv;gv(7:9)=vv;

	x1=0.5853132347709715;y1=0.2073433826145142;z1=1-x1-y1;  % multiplicity 3
 	w1=0.2164790926342230/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(10:12)=xv;gcy(10:12)=yv;gv(10:12)=vv;

	x1=0.1221843885990187;y1=0.4389078057004907;z1=1-x1-y1;  % multiplicity 3
 	w1=0.20798741611661160/3;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,3);
	gcx(13:15)=xv;gcy(13:15)=yv;gv(13:15)=vv;

	x1=0;y1=0.858870281282263640;z1=1-x1-y1;  % multiplicity 6
	w1=0.04417430269980344/6;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,6);
	gcx(16:21)=xv;gcy(16:21)=yv;gv(16:21)=vv;

	x1=0.04484167758913055;y1=0.6779376548825902;z1=1-x1-y1;  % multiplicity 6
	w1=0.2463378925757316/6;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,6);
	gcx(22:27)=xv;gcy(22:27)=yv;gv(22:27)=vv;

	x1=1/3;y1=1/3;z1=1-x1-y1;  % multiplicity 1
	w1=0.08797730116222190;
	[xv,yv,vv]=multiplicity(x1,y1,w1,xc,yc,1);
	gcx(end)=xv;gcy(end)=yv;gv(end)=vv;

otherwise
	disp('Polynomial degrees are 1,2,3,4,5,7,8,11 in trigauc')
	gcx=[];gcy=[];gv=[];
end


function [xv,yv,vv]=multiplicity(x,y,v,xc,yc,m)
z=[x,y,1-x-y];
switch(m)
case 1
	xv=[xc(1)*z(1)+xc(2)*z(2)+xc(3)*z(3)];
	yv=[yc(1)*z(1)+yc(2)*z(2)+yc(3)*z(3)];
    vv=v;
case 3
	xv=[xc(1)*z(1)+xc(2)*z(2)+xc(3)*z(3),xc(1)*z(2)+xc(2)*z(3)+xc(3)*z(1),...
	xc(1)*z(3)+xc(2)*z(1)+xc(3)*z(2)];
	yv=[yc(1)*z(1)+yc(2)*z(2)+yc(3)*z(3),yc(1)*z(2)+yc(2)*z(3)+yc(3)*z(1),...
	yc(1)*z(3)+yc(2)*z(1)+yc(3)*z(2)];
	vv=v*ones(1,3);
case 6
	xv=[xc(1)*z(1)+xc(2)*z(2)+xc(3)*z(3),xc(1)*z(1)+xc(2)*z(3)+xc(3)*z(2),...
	xc(1)*z(2)+xc(2)*z(1)+xc(3)*z(3),xc(1)*z(2)+xc(2)*z(3)+xc(3)*z(1),...
	xc(1)*z(3)+xc(2)*z(1)+xc(3)*z(2),xc(1)*z(3)+xc(2)*z(2)+xc(3)*z(1)];
	yv=[yc(1)*z(1)+yc(2)*z(2)+yc(3)*z(3),yc(1)*z(1)+yc(2)*z(3)+yc(3)*z(2),...
	yc(1)*z(2)+yc(2)*z(1)+yc(3)*z(3),yc(1)*z(2)+yc(2)*z(3)+yc(3)*z(1),...
	yc(1)*z(3)+yc(2)*z(1)+yc(3)*z(2),yc(1)*z(3)+yc(2)*z(2)+yc(3)*z(1)];
	vv=v*ones(1,6);
end
