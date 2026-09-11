function results = verify_korn_constants(outdir)
% Verify the Korn-constant values quoted in Remark A.6 on the box
% (0,3)x(0,1)x(0,1) with GammaD = {y=0} u {y=1} u {x=0}:
%   broken     max ||grad_h v||^2 / Xi(v)^2 over discontinuous P1 fields,
%   conforming max ||grad u||^2 / (||eps^dev u||^2 + ||u||^2_GammaD) over P1,
%   H^1_0      max ||grad u||^2 / ||eps^dev u||^2 over P1 vanishing on the
%              whole boundary (exact value 2, approached from below).
% The broken problem is factored once with a nested dissection Cholesky
% factor; the 36x12x12 case needs about 6 GB of memory.

if(nargin<1 || isempty(outdir)), outdir=fullfile(pwd,'results'); end
if(~exist(outdir,'dir')), mkdir(outdir); end

m=[9;18;27;36];n=numel(m);
[broken,edge_energy_share,conforming,h10]=deal(nan(n,1));
for k=1:n
    [A,B,tet,x,~,z]=broken_forms(m(k),m(k)/3,m(k)/3);
    nd=size(A,1);q=dissect(A);[Lf,flag]=chol(A(q,q),'lower');assert(flag==0)
    Bq=B(q,q);fun=@(v) Lf\(Bq*(Lf'\v));
    [V,D]=eigs(fun,nd,1,'largestabs','IsFunctionSymmetric',true, ...
        'Tolerance',1e-10,'MaxIterations',5000);
    broken(k)=D(1,1);
    xv=zeros(nd,1);xv(q)=Lf'\V(:,1);nele=size(tet,1);
    eT=sum(reshape(xv.*(B*xv),12,nele),1)';eT=eT/sum(eT);
    xc=mean(x(tet),2);zc=mean(z(tet),2);h=3/m(k);
    edge=xc>3-2*h&(zc<2*h|zc>1-2*h);edge_energy_share(k)=sum(eT(edge));
    clear A B Lf Bq V fun
    [conforming(k),h10(k)]=conforming_constants(m(k),m(k)/3,m(k)/3,k<n);
end
mesh=m;
results=table(mesh,broken,edge_energy_share,conforming,h10);
writetable(results,fullfile(outdir,'korn_constants.csv'));

check_close('broken Korn constants',broken,[6.3742;11.1434;16.0720;20.7317],2e-3,0);
check_close('edge energy shares',edge_energy_share(2:2:4),[0.809;0.647],2e-3,0);
check_close('conforming Korn constants',conforming,[6.818;10.870;14.366;17.209],2e-3,0);
check_close('H10 quotients',h10(1:3),[1.9888;1.9976;1.9990],5e-4,0);
fprintf('Verified the Korn constants of Remark A.6 in %s\n',outdir)
end

function [cc,ch]=conforming_constants(nx,ny,nz,do_h10)
tol=1e-9;Dtil=blkdiag(eye(3)-ones(3)/3,0.5*eye(3));
[tet,x,y,z]=meshbox(3,1,1,nx,ny,nz);nele=size(tet,1);nd=3*numel(x);
ra=zeros(nele*144,1);ca=ra;va=ra;vb=ra;u=0;
for iel=1:nele
    iv=tet(iel,:);[~,fx,fy,fz,V]=basis3(mean(x(iv)),mean(y(iv)),mean(z(iv)),x(iv),y(iv),z(iv));
    Bm=strain_matrix(fx,fy,fz);G=[fx fy fz];
    Ae=V*(Bm'*Dtil*Bm);Be=kron(V*(G*G'),eye(3));
    eqs=reshape(3*(iv-1)+(1:3)',[],1);[I,J]=ndgrid(eqs,eqs);
    ra(u+1:u+144)=I(:);ca(u+1:u+144)=J(:);va(u+1:u+144)=Ae(:);vb(u+1:u+144)=Be(:);u=u+144;
end
A=sparse(ra,ca,va,nd,nd);B=sparse(ra,ca,vb,nd,nd);A=(A+A')/2;B=(B+B')/2;
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');cnt=accumarray(ic,1);
v1=[x(fu(:,2))-x(fu(:,1)),y(fu(:,2))-y(fu(:,1)),z(fu(:,2))-z(fu(:,1))];
v2=[x(fu(:,3))-x(fu(:,1)),y(fu(:,3))-y(fu(:,1)),z(fu(:,3))-z(fu(:,1))];
far=0.5*sqrt(sum(cross(v1,v2,2).^2,2));X=x(fu);Y=y(fu);
isD=find((cnt==1)&(all(abs(Y)<tol,2)|all(abs(Y-1)<tol,2)|all(abs(X)<tol,2)));
rm=zeros(numel(isD)*81,1);cm=rm;vm=rm;um=0;M0=(eye(3)+ones(3))/12;
for e=isD'
    for c=1:3
        dofs=3*(fu(e,:)'-1)+c;[I,J]=ndgrid(dofs,dofs);
        rm(um+1:um+9)=I(:);cm(um+1:um+9)=J(:);vm(um+1:um+9)=far(e)*M0(:);um=um+9;
    end
end
Mb=sparse(rm(1:um),cm(1:um),vm(1:um),nd,nd);Ac=(A+Mb+(A+Mb)')/2;
cc=eigs(B,Ac,1,'largestabs','Tolerance',1e-8,'MaxIterations',3000);
ch=NaN;
if(do_h10)
    onB=abs(x)<tol|abs(x-3)<tol|abs(y)<tol|abs(y-1)<tol|abs(z)<tol|abs(z-1)<tol;
    keep=true(nd,1);for c=1:3, keep(3*(find(onB)-1)+c)=false; end
    ch=eigs(B(keep,keep),A(keep,keep),1,'largestabs','Tolerance',1e-8,'MaxIterations',3000);
end
end

function [A,B,tet,x,y,z]=broken_forms(nx,ny,nz)
% A: sum_T ||eps^dev(v)||_T^2 + sum_F h_F^-1 ||[v]||_F^2 (interior and GammaD)
% B: sum_T ||grad v||_T^2 ; v discontinuous P1, dof 12(T-1)+3(i-1)+c
tol=1e-9;Dtil=blkdiag(eye(3)-ones(3)/3,0.5*eye(3));
[tet,x,y,z]=meshbox(3,1,1,nx,ny,nz);nele=size(tet,1);nd=12*nele;
ra=zeros(nele*144,1);ca=ra;va=ra;vb=ra;u=0;vol=zeros(nele,1);
for iel=1:nele
    iv=tet(iel,:);[~,fx,fy,fz,V]=basis3(mean(x(iv)),mean(y(iv)),mean(z(iv)),x(iv),y(iv),z(iv));
    vol(iel)=V;Bm=strain_matrix(fx,fy,fz);G=[fx fy fz];
    Ae=V*(Bm'*Dtil*Bm);Be=kron(V*(G*G'),eye(3));
    eqs=12*(iel-1)+(1:12)';[I,J]=ndgrid(eqs,eqs);
    ra(u+1:u+144)=I(:);ca(u+1:u+144)=J(:);va(u+1:u+144)=Ae(:);vb(u+1:u+144)=Be(:);u=u+144;
end
A=sparse(ra,ca,va,nd,nd);B=sparse(ra,ca,vb,nd,nd);clear ra ca va vb
F=[tet(:,[2,3,4]);tet(:,[1,3,4]);tet(:,[1,2,4]);tet(:,[1,2,3])];
[fu,~,ic]=unique(sort(F,2),'rows');clear F
nface=size(fu,1);cnt=accumarray(ic,1);
kel=mod((0:4*nele-1)',nele)+1;Tk=tet(kel,:);Fk=fu(ic,:);loc=zeros(4*nele,3);
for mm=1:3, [~,loc(:,mm)]=max(Tk==Fk(:,mm),[],2); end
clear Tk Fk
[ics,ord]=sort(ic);starts=[1;find(diff(ics))+1];
first=ord(starts);second=zeros(nface,1);second(cnt==2)=ord(starts(cnt==2)+1);
v1=[x(fu(:,2))-x(fu(:,1)),y(fu(:,2))-y(fu(:,1)),z(fu(:,2))-z(fu(:,1))];
v2=[x(fu(:,3))-x(fu(:,1)),y(fu(:,3))-y(fu(:,1)),z(fu(:,3))-z(fu(:,1))];
far=0.5*sqrt(sum(cross(v1,v2,2).^2,2));
hF=accumarray(ic,vol(kel),[nface,1])./(2*far);X=x(fu);Y=y(fu);
isD=(cnt==1)&(all(abs(Y)<tol,2)|all(abs(Y-1)<tol,2)|all(abs(X)<tol,2));
M0=(eye(3)+ones(3))/12;sgn=[1;-1];
rj=zeros(nface*108,1);cj=rj;vj=rj;uj=0;
for e=find((cnt==2)|isD)'
    if cnt(e)==2, occ=[first(e);second(e)]; else, occ=first(e); end
    ns=numel(occ);sg=kron(sgn(1:ns),ones(3,1));
    blk=(far(e)/hF(e))*(sg*sg').*repmat(M0,ns,ns);
    pos=12*(kel(occ)-1)+3*(loc(occ,:)-1);
    for c=1:3
        dofs=reshape((pos+c)',[],1);[I,J]=ndgrid(dofs,dofs);n2=numel(I);
        rj(uj+1:uj+n2)=I(:);cj(uj+1:uj+n2)=J(:);vj(uj+1:uj+n2)=blk(:);uj=uj+n2;
    end
end
A=A+sparse(rj(1:uj),cj(1:uj),vj(1:uj),nd,nd);A=(A+A')/2;B=(B+B')/2;
end

function Bm=strain_matrix(fx,fy,fz)
Bm=zeros(6,12);
for i=1:4
    cx=3*(i-1);
    Bm(1,cx+1)=fx(i);Bm(2,cx+2)=fy(i);Bm(3,cx+3)=fz(i);
    Bm(4,cx+1)=fy(i);Bm(4,cx+2)=fx(i);
    Bm(5,cx+1)=fz(i);Bm(5,cx+3)=fx(i);
    Bm(6,cx+2)=fz(i);Bm(6,cx+3)=fy(i);
end
end

function check_close(name,actual,reference,atol,rtol)
err=abs(actual-reference);lim=atol+rtol*abs(reference);
if(any(err>lim))
    error('reproduce:regression','%s failed: max error %.3e',name,max(err))
end
end
