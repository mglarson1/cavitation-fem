function results=verify_stokes_qp_reference(outdir,meshes)
% Independent unsmoothed reference and KKT/working-set diagnostics.
if nargin<1,outdir=fullfile(pwd,'results');end
if nargin<2,meshes=[128,8;256,16];end
if ~exist(outdir,'dir'),mkdir(outdir);end
rows=[];
for k=1:size(meshes,1)
 nx=meshes(k,1);nz=meshes(k,2);[tri,x,y,di,dv,Hf]=pit_mesh(nx,nz);
 for tol=[1e-10,1e-12]
 for ktol=[1e-10,1e-12,1e-14]
  tic;[u,p,info]=solve_stokes_qp_reference(tri,x,y,di,dv,tol,ktol);elapsed=toc;
  xs=linspace(0,24,8001)';TR=triangulation(tri,x,y);ei=pointLocation(TR,[xs,1-0.5*Hf(xs)]);w=xs>4&xs<20;ei=ei(w);
  peak=max(info.qp_p(ei));pd=info.polished_div;pp=info.polished_p;
  cavity=0.003*nnz(~info.polished_act(ei));
  rows=[rows;nx,nz,size(tri,1),tol,ktol,info.polished_iterations,info.exitflag,info.iterations,elapsed,peak, ...
   info.stationarity,info.min_pressure,info.min_divergence,info.complementarity, ...
   max(pp(ei)),cavity,info.polished_stationarity,min(pp),min(pd), ...
   info.polished_complementarity,info.polished_energy-info.energy]; %#ok<AGROW>
  results=array2table(rows,'VariableNames',{'nx','nz','elements','tolerance', ...
   'kkt_tolerance','polish_iterations','exitflag','iterations','seconds','peak_qp','stationarity_qp','min_pressure_qp', ...
   'min_divergence_qp','complementarity_qp','peak_polished','cavity_polished', ...
   'stationarity_polished','min_pressure_polished','min_divergence_polished', ...
   'complementarity_polished','energy_change'});
  writetable(results,fullfile(outdir,'stokes_qp_reference.csv'));
  save(fullfile(outdir,sprintf('stokes_qp_%d_%g_kkt_%g.mat',nx,tol,ktol)),'u','p','info','tri','x','y','di','dv');
  fprintf('QP %d elements tol=%g: peak=%.9g cavity=%.6f; polished min p=%.3e div=%.3e stat=%.3e, time=%.1fs\n', ...
   size(tri,1),tol,peak,cavity,min(pp),min(pd),info.polished_stationarity,elapsed);
 end
 end
end
end

function [tri,x,z,dind,dval,Hf]=pit_mesh(nx,nz)
Hf=@(x) 1+exp(-(x-12).^2);
[tri,x,zh]=meshrect(24,1,nx,nz);z=1-zh.*Hf(x);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
a=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
b=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*a-1;2*a;2*b-1;2*b];dval=[ones(numel(a),1);zeros(numel(a)+2*numel(b),1)];
end
