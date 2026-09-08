function results=verify_pressure_end_flat_channel(outdir)
% Exact Couette checks: tangential velocity zero, prescribed normal traction.
if nargin<1,outdir=fullfile(pwd,'results','pressure_ends_2026-09-08');end
if ~exist(outdir,'dir'),mkdir(outdir);end
rows=[];
for nz=[8,16]
 nx=16*nz;[tri,x,z]=meshrect(24,1,nx,nz);
 E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];[ed,~,ic]=unique(sort(E,2),'rows');bd=accumarray(ic,1)==1;
 slide=find(bd&z(ed(:,1))==1&z(ed(:,2))==1);wall=find(bd&z(ed(:,1))==0&z(ed(:,2))==0);
 le=find(bd&x(ed(:,1))==0&x(ed(:,2))==0);ri=find(bd&x(ed(:,1))==24&x(ed(:,2))==24);ends=[le;ri];
 di=[2*slide-1;2*slide;2*wall-1;2*wall;2*ends];dv=[ones(numel(slide),1);zeros(numel(slide)+2*numel(wall)+numel(ends),1)];
 dc=false(size(ed,1),2);dc([slide;wall],:)=true;dc(ends,2)=true;
 uexact=zeros(2*size(ed,1),1);uexact(1:2:end)=mean(z(ed),2);
 for pb=[0,1]
  force=zeros(size(uexact));force(2*le-1)=pb*abs(diff(z(ed(le,:)),1,2));force(2*ri-1)=-pb*abs(diff(z(ed(ri,:)),1,2));
  op=struct('gamma1',1,'lamfac',-2/3,'solver','qp','verbose',0, ...
   'dirichlet_components',dc,'load',force,'reference_u',uexact, ...
   'reference_p',pb*ones(size(tri,1),1));
  [ux,uz,p,info]=solvedisccrs(tri,x,z,100,di,dv,op);
  eu=max(abs(ux-uexact(1:2:end)));ev=max(abs(uz));ep=max(abs(p-pb));
  assert(max([eu,ev,ep,info.reference_stationarity,info.reference_divergence])<1e-10);
  rows=[rows;nx,nz,pb,eu,ev,ep,info.stationarity,info.reference_stationarity, ...
    info.reference_divergence,info.qp.qp_iterations,info.iterations]; %#ok<AGROW>
  fprintf('Flat nx%d nz%d p_ext%g: errors ux%.3e uz%.3e p%.3e; exact residual%.3e\n',nx,nz,pb,eu,ev,ep,info.reference_stationarity);
 end
end
results=array2table(rows,'VariableNames',{'nx','nz','external_pressure','ux_error', ...
 'uz_error','pressure_error','stationarity','exact_stationarity','exact_divergence','qp_iterations','polish_iterations'});
writetable(results,fullfile(outdir,'flat_channel.csv'));
% Reynolds zero-pressure flat channel: no forcing and homogeneous end data.
[tri,x,y]=meshrect(24,4*24/256,256,4);dn=find(x==0|x==24);
[P,lam,ir]=solvereynolds(tri,x,y,@(x,y)1,@(x,y)0,1,dn);
assert(ir.converged && max(abs(P))<1e-13 && max(abs(lam))<1e-13);
end
