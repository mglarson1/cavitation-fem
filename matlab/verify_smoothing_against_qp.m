function results=verify_smoothing_against_qp(outdir,svals,tolerances,filename)
% Compare the actual smoothed fields and cavity diagnostic with a KKT-checked QP.
if nargin<1,outdir=fullfile(pwd,'results');end
if nargin<2,svals=[1e-4,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16,1e-18,1e-20];end
if nargin<3,tolerances=1e-14;end
if nargin<4,filename='smoothing_against_qp.csv';end
ref=load(fullfile(outdir,'stokes_qp_256_1e-12_kkt_1e-12.mat'));
tri=ref.tri;x=ref.x;y=ref.y;area=ref.info.area;B=ref.info.B;K=ref.info.K;
xs=linspace(0,24,8001)';ei=pointLocation(triangulation(tri,x,y), ...
 [xs,1-.5*(1+exp(-(xs-12).^2))]);ei=ei(xs>4&xs<20);
c0=.003*nnz(~ref.info.polished_act(ei));p0=ref.p;u0=ref.u;
rows=[];
for s=svals
 for tol=tolerances
 opts=struct('maxit',40,'tol',tol,'verbose',0,'damp',0);
 [ux,uy,p,info]=solvedisccrn(tri,x,y,100,s,ref.di,ref.dv,opts);
 u=zeros(size(u0));u(1:2:end)=ux;u(2:2:end)=uy;du=u-u0;
 cavity=.003*nnz(p(ei)<sqrt(s/100));D=B*u./area;
 station=norm(K(ref.info.free,:)*u-B(:,ref.info.free)'*p,inf);
 rows=[rows;s,tol,info.converged,info.iterations,info.relres,max(p(ei)),cavity,c0, ...
  abs(cavity-c0)/c0,sqrt(sum(area.*(p-p0).^2)/sum(area.*p0.^2)), ...
  sqrt(max(0,du'*K*du)/(u0'*K*u0)),min(p),min(D), ...
  max(abs(p.*D-s/100)),station,max(abs(D-info.div)),.003*nnz(xor(p(ei)<sqrt(s/100),~ref.info.polished_act(ei)))]; %#ok<AGROW>
 results=array2table(rows,'VariableNames',{'smoothing','tolerance','converged','iterations', ...
  'relative_residual','peak','cavity','reference_cavity','relative_cavity_error', ...
  'relative_pressure_l2_error','relative_velocity_energy_error','min_pressure', ...
  'min_divergence','central_path_error','stationarity','divergence_assembly_difference','classification_mismatch_length'});
 writetable(results,fullfile(outdir,filename));
 fprintf('s=%g converged%d its%d cavity%.6f error%.2f%% p_err%.3e minD%.3e\n', ...
  s,info.converged,info.iterations,cavity,100*rows(end,9),rows(end,10),min(D));
 end
end
end
