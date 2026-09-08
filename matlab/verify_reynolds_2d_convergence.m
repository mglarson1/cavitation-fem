function results=verify_reynolds_2d_convergence(outdir,ns)
% Raw L2 and H1 errors on uniformly refined, shape-regular square meshes.
% Exact obstacle solution P=(1/4-x^2-y^2)_+^2; d=1 on (-1,1)^2.
% f=2-16r^2 inside r<1/2 and f=-2 outside. lambda=2 outside.
if nargin<1,outdir=fullfile(pwd,'results');end
if nargin<2,ns=[16;32;64;128];end
if ~exist(outdir,'dir'),mkdir(outdir);end
rows=[];
for n=ns(:)'
 [tri,x,y]=meshrect(2,2,n,n);x=x-1;y=y-1;
 df=@(x,y) 1+0*x+0*y;
 ff=@(x,y) max(2-16*(x.^2+y.^2),-2);
 tic;evalc('[P,lam,info]=solvereynolds(tri,x,y,df,ff,100);');
 assert(info.converged&&info.signok)
 [e2,e1]=raw_errors(tri,x,y,P,6);
 [e2check,e1check]=raw_errors(tri,x,y,P,10);
 free=~info.isdir;kkt=norm(info.K(free,:)*P-info.F(free)-info.m(free).*lam(free),inf);
 rows=[rows;n,size(tri,1),numel(x),2*sqrt(2)/n,e2check,e1check, ...
  abs(e2-e2check)/e2check,abs(e1-e1check)/e1check,info.iterations, ...
  min(P),min(lam),info.complementarity,kkt,toc]; %#ok<AGROW>
 results=array2table(rows,'VariableNames',{'n','elements','nodes','hmax','l2_error', ...
  'h1_seminorm_error','l2_quadrature_relative_change','h1_quadrature_relative_change', ...
  'iterations','min_pressure','min_multiplier','complementarity','equilibrium','seconds'});
 results.l2_rate=[NaN;log(results.l2_error(1:end-1)./results.l2_error(2:end))/log(2)];
 results.h1_rate=[NaN;log(results.h1_seminorm_error(1:end-1)./results.h1_seminorm_error(2:end))/log(2)];
 writetable(results,fullfile(outdir,'reynolds_2d_convergence.csv'));
 writetable(results,fullfile(outdir,'table02_reynolds_convergence.csv'));
 fprintf('Reynolds n=%d, L2=%.8e H1=%.8e rates %.4f %.4f quadrature changes %.2e %.2e\n', ...
  n,e2check,e1check,results.l2_rate(end),results.h1_rate(end),rows(end,7),rows(end,8));
end
end

function [l2,h1]=raw_errors(tri,x,y,P,nq)
% Composite Duffy-Gauss integration; cut triangles receive two refinement levels (16 subtriangles).
% The raw affine FE polynomial is evaluated on every integration subtriangle.
[q,w]=gauss01(nq);[a,b]=ndgrid(q,q);[wa,wb]=ndgrid(w,w);
l=[1-a(:),a(:).*(1-b(:)),a(:).*b(:)];weights=2*wa(:).*wb(:).*a(:);
nt=size(tri,1);xv=x(tri);yv=y(tri);pv=P(tri);
detj=(xv(:,2)-xv(:,1)).*(yv(:,3)-yv(:,1))-(xv(:,3)-xv(:,1)).*(yv(:,2)-yv(:,1));
gx=sum(pv.*[yv(:,2)-yv(:,3),yv(:,3)-yv(:,1),yv(:,1)-yv(:,2)],2)./detj;
gy=sum(pv.*[xv(:,3)-xv(:,2),xv(:,1)-xv(:,3),xv(:,2)-xv(:,1)],2)./detj;
% Conservative circle-intersection classification using a bounding box.
rmin=max(0,max(min(xv,[],2),-max(xv,[],2))).^2 + ...
 max(0,max(min(yv,[],2),-max(yv,[],2))).^2;
rmax=max(xv.^2+yv.^2,[],2);cut=rmin<.25&rmax>.25;
le=0;he=0;
for t=1:nt
 xx=xv(t,:);yy=yv(t,:);
 if cut(t)
  % Two uniform subdivisions (16 subtriangles) control the nonpolynomial cut.
  for depth=1:2
   x12=(xx(:,1)+xx(:,2))/2;x23=(xx(:,2)+xx(:,3))/2;x31=(xx(:,3)+xx(:,1))/2;
   y12=(yy(:,1)+yy(:,2))/2;y23=(yy(:,2)+yy(:,3))/2;y31=(yy(:,3)+yy(:,1))/2;
   xx=[xx(:,1),x12,x31;x12,xx(:,2),x23;x31,x23,xx(:,3);x12,x23,x31];
   yy=[yy(:,1),y12,y31;y12,yy(:,2),y23;y31,y23,yy(:,3);y12,y23,y31];
  end
 end
 area=abs((xx(:,2)-xx(:,1)).*(yy(:,3)-yy(:,1))-(xx(:,3)-xx(:,1)).*(yy(:,2)-yy(:,1)))/2;
 X=l*xx';Y=l*yy';qval=max(.25-X.^2-Y.^2,0);
 exact=qval.^2;ex=-4*X.*qval;ey=-4*Y.*qval;
 uh=pv(t,1)+gx(t)*(X-xv(t,1))+gy(t)*(Y-yv(t,1));
 le=le+weights'*(exact-uh).^2*area;
 he=he+weights'*((ex-gx(t)).^2+(ey-gy(t)).^2)*area;
end
l2=sqrt(le);h1=sqrt(he);
end

function [x,w]=gauss01(n)
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);[V,D]=eig(diag(b,1)+diag(b,-1));
[x,i]=sort(diag(D));x=(x+1)/2;w=(V(1,i).^2)';
end
