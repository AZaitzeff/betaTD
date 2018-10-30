function [mapall,dict,kappa]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,dx,dy,DT,dtstop,numsub)
if nargin<8
    dx=1/100;
    dy=1/100;
end

if nargin<9
    DT=.02;
end
if nargin<10
    dtstop=2^(-12);
end
if nargin<11
    numsub=200;
end
dt=DT;

T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
K=max(mapall(:));
[m,n]=size(mapall);%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[m*n,1]);
u=cell(1,K);
S=cell(1,K);
MAXITER=1000;
for k=1:K
    u{k}=mapall==k;
    S{k}=1-u{k}*2;
end

changeu=u;
lastu=u;
for t=1:MAXITER

[ls]=td2dz(u,dt,dx,dy);
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

if totalnum<2
    if dt<dtstop
        break
    else
        dt=dt/2;
    end
end
lastu=u;


end

for k=1:K
    mapall(u{k}>0)=k;
end
