function combined=audit_domain_length_results(outdir)
% Verify the control variables and accepted nonlinear solutions from MAT files.
if nargin<1,outdir=fullfile(pwd,'results','domain_length_2026-09-08');end
files={'domain_0_24.mat','one_sided/domain_-12_24.mat', ...
 'one_sided/domain_0_36.mat','domain_-12_36.mat','domain_-36_60.mat'};
refS=[];refR=[];
for k=1:numel(files)
 q=load(fullfile(outdir,files{k}));
 a=q.x(q.tri);b=q.z(q.tri);sel=min(a,[],2)>=0 & max(a,[],2)<=24;
 geomS=sortrows([a(sel,:),b(sel,:)]);
 a=q.xR(q.triR);b=q.yR(q.triR);sel=min(a,[],2)>=0 & max(a,[],2)<=24;
 geomR=sortrows([a(sel,:),b(sel,:)]);
 if isempty(refS),refS=geomS;refR=geomR;end
 assert(isequal(refS,geomS),'The central Stokes triangles changed');
 assert(isequal(refR,geomR),'The central Reynolds triangles changed');
 assert(q.iR.converged && q.iS.converged && q.iR.signok && q.iS.signok);
 assert(q.iR.complementarity<1e-18 && q.iS.compl<5e-13);
 fprintf('%s: %d Stokes and %d Reynolds central triangles identical; nonlinear checks passed.\n', ...
  files{k},size(geomS,1),size(geomR,1));
end
sym=readtable(fullfile(outdir,'domain_length_sensitivity.csv'));
side=readtable(fullfile(outdir,'one_sided','domain_length_sensitivity.csv'));
assert(isequal(sym{1,13:26},side{1,13:26}),'Repeated baseline changed');
combined=[sym(1,:);side(2:3,:);sym(2:3,:)];
writetable(combined,fullfile(outdir,'table12_domain_length.csv'));
end
