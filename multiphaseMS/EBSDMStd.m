function [u,dict,kappa]=EBSDMStd(u,EBSD,CI,dict,kappa,fid,DT)
if nargin<7
    DT=.02;
end
dt=DT;
numsub=200;
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
K=size(u,2);
[m,n]=size(u{1});%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[m*n,1]);
S=u;
changeu=u;
lastu=u;
MAXITER=1000;
for k=1:K
    S{k}(:)=CIflat.*alpbmetric(EBSDflat,dict{k})';
end
for t=1:MAXITER


[ls]=td2dz(u,dt);
phi=zeros(m,n,K);
for k=1:K
    phi(:,:,k)=2/(sqrt(dt))*(1-ls{k})+fid*S{k};
end
for k=1:K
    temp=phi;
    temp(:,:,k)=[];
    tempmin=squeeze(max(-temp,[],3));
    u{k}=(squeeze(-phi(:,:,k))-tempmin);
    
end
totalnum=0;
for k=1:K
    ind=(u{k}(:)>0);
    totalnum=totalnum+sum(abs((lastu{k}(:)>0)-ind));
    if sum(ind)>10
    num=sum(abs((changeu{k}(:)>0)-ind))/sum(ind);
    
    if (num)>.05
        changeu{k}=u{k};
        mask1=u{k}>0;
        indices=find(mask1(:));
        if sum(CIflat(indices))>1e-4
            if numel(indices)>numsub
                newind=datasample(indices,numsub,'Weights',CIflat(indices));
                EBSDtemp=EBSDflat(newind,:);
                CItemp=ones(size(CIflat(newind)));
            else
                EBSDtemp=EBSDflat(indices,:);
                CItemp=CIflat(indices);
            end
            [newg1, ~, ~, ~] = VMFEM(EBSDtemp, Pall,CItemp,1,3,dict{k},kappa{k});
            dict{k}=newg1;
            S{k}(:)=CIflat.*alpbmetric(EBSDflat,dict{k})';
        end
    end
    end
end
if totalnum<1
    dt=dt/2;
    if dt<1/2^12
        break
    end
end
lastu=u;


end
