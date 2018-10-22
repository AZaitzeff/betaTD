function [u,dict,kappa]=EBSDMStd(u,EBSD,CI,dict,kappa,fid,DT,numsub)
if nargin<7
    DT=.02;
end
if nargin<8
    numsub=120;
end
dt=DT;

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
S=cell(1,K);
changeu=u;
lastu=u;
MAXITER=1000;
for k=1:K
    S{k}=1-u{k}*2;
end
for t=1:MAXITER

[ls]=td2dz(u,dt);
phi=zeros(m,n,K);
for k=1:K
    mask=(ls{k}<.99)&(ls{k}>.01)&((S{k}>.95)|(S{k}<-.01));
    S{k}(mask)=CIflat(mask).*alpbmetric(EBSDflat(mask,:),dict{k})';
end

for k=1:K
    phi(:,:,k)=2/(sqrt(dt))*(1-ls{k})+fid*S{k};
end

[~,argmin]=min(phi,[],3);
for k=1:K
    u{k}=argmin==k;
    
end
totalnum=0;
k=1;
while k<=K
    ind=(u{k}(:)>0);
    totalnum=totalnum+sum(abs((lastu{k}(:)>0)-ind));
    if sum(ind)>5
    num=sum(abs((changeu{k}(:)>0)-ind))/sum(ind);
    
    if (num)>=.4
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
            [newg1, kap, ~, ~] = VMFEM(EBSDtemp, Pall,CItemp,1,1,dict{k},kappa{k});
            dict{k}=newg1;
            kappa{k}=kap;
            S{k}=1-u{k}*2;
        end
    end
    k=k+1;
    else
        dict{k}=dict{K};
        u{k}=u{K};
        kappa{k}=kappa{K};
        kappa(K)=[];
        u(K)=[];
        dict(K)=[];
        changeu{k}=changeu{K};
        lastu{k}=lastu{K};
        S{k}=S{K};
        S(K)=[];
        K=K-1;
    end
end

if totalnum<100
    if dt<1/2^12
        if totalnum<2
            break
        end
    else
        dt=dt/2;
    end
end
lastu=u;


end
