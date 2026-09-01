function [tri,xnod,ynod,bnod]=sqtriKb(tril,xnodl,ynodl,edges,bnodl)

%       6
%      | \
%      |  \
%      |   \
%      |    \
%      4     5
%      |      \
%      |       \
%      |        \
%     1 --- 2 -- 3
%
nno=length(xnodl);
nedge=size(edges,1);
xnod=zeros(nno+nedge,1);ynod=xnod;bnod=xnod;
ynod(1:nno)=ynodl;xnod(1:nno)=xnodl;bnod(1:nno)=bnodl;
nele=size(tril,1);
tri=zeros(nele,6);ed=zeros(3,1);

nn=[1,2;2,3;3,1];
for iel=1:nele
    ivv=tril(iel,:);
    for j=1:3
        ind=nn(j,:);k=sort(ivv(ind));
        ed(j)=find(edges(:,1)==k(1) & edges(:,2)==k(2));
        xnod(nno+ed(j))=(xnod(k(1))+xnod(k(2)))/2;
        ynod(nno+ed(j))=(ynod(k(1))+ynod(k(2)))/2;
        if(bnod(k(1))==1 & bnod(k(2))==1)
            bnod(nno+ed(j))=1;
        end
    end
    ivl=[ivv(1),nno+ed(1),ivv(2),nno+ed(3),nno+ed(2),ivv(3)];
    tri(iel,:)=ivl;
end


