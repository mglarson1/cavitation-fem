function results = verify_mechanical_pressure(outdir)
% Verify the mechanical-pressure numbers quoted in Section 5.6: steepest pit,
% stabilized Crouzeix--Raviart, traction-free ends, mesh 256x16.  With first
% Lame parameter lamfac*mu the mechanical pressure is
%   pbar = -tr(sigma)/3 = p - (2/3 + lamfac) mu div u ,
% so pbar = p for the deviatoric law (lamfac = -2/3) and
% pbar = p - (2 mu/3) div u for the customary law (lamfac = 0).

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

L=24;xp=12;c=1;V=1;nx=256;nz=16;
Hf=@(x) c*(1+exp(-(x-xp).^2));
[tri,xh,zh]=meshrect(L,1,nx,nz);x=xh;z=c-zh.*Hf(xh);
E=[tri(:,[2,3]);tri(:,[1,3]);tri(:,[1,2])];
[eu,~,ic]=unique(sort(E,2),'rows');onb=(accumarray(ic,1)==1);
sl=find(onb&abs(zh(eu(:,1)))<1e-9&abs(zh(eu(:,2)))<1e-9);
sh=find(onb&abs(zh(eu(:,1))-1)<1e-9&abs(zh(eu(:,2))-1)<1e-9);
dind=[2*sl-1;2*sl;2*sh-1;2*sh];
dval=[V*ones(numel(sl),1);zeros(numel(sl),1);zeros(2*numel(sh),1)];
TR=triangulation(tri,x,z);xs=linspace(0,L,8001)';
eid=pointLocation(TR,[xs,c-0.5*Hf(xs)]);ok=~isnan(eid);
w=xs>4*c&xs<L-4*c;dx=xs(2)-xs(1);

lamfac=[-2/3;0];law=["deviatoric";"customary"];
[cavitated_elements,negative_fraction,min_pbar_cavity,cavity_length,iterations]=deal(zeros(2,1));
for k=1:2
    op=struct('gamma1',1,'lamfac',lamfac(k),'verbose',0,'solver','as');
    [~,~,p,info]=solvedisccrs(tri,x,z,100,dind,dval,op);
    assert(info.converged && info.signok)
    cav=~info.act;pbar=p-(2/3+lamfac(k))*info.div;
    cs=false(size(xs));cs(ok)=cav(eid(ok));
    cavitated_elements(k)=nnz(cav);iterations(k)=info.iterations;
    negative_fraction(k)=nnz(pbar(cav)<-1e-12)/nnz(cav);
    min_pbar_cavity(k)=min(pbar(cav));cavity_length(k)=dx*nnz(cs&w);
end
results=table(law,lamfac,iterations,cavitated_elements,negative_fraction, ...
    min_pbar_cavity,cavity_length);
writetable(results,fullfile(outdir,'prose_mechanical_pressure.csv'));

assert(cavitated_elements(1)==492 && cavitated_elements(2)==676)
assert(negative_fraction(1)==0 && min_pbar_cavity(1)>=-1e-12)
assert(negative_fraction(2)==1 && abs(min_pbar_cavity(2)+0.64)<5e-3)
assert(abs(cavity_length(1)-0.843)<2e-4 && abs(cavity_length(2)-1.314)<2e-4)
fprintf('Verified the mechanical-pressure claims in %s\n',outdir)
end
