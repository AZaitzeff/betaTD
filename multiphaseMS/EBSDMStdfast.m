function [mapall,dict,kappa]=EBSDMStdfast(mapall,EBSD,CI,dict,kappa,fid,dx,dy,DT,rmspots,mexed,numsub)
if nargin<8
    dx=1/100;
    dy=1/100;
end
if nargin<9
    DT=.02;
end
if nargin<12
    numsub=200;
end
if nargin<11
    mexed=0;
end
if nargin<10
    rmspots=0;
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
MAXITER=1000;
changenum=zeros(1,K);
curmin=ones(m,n)*fid*2;
active=ones(1,K);
for t=1:MAXITER
[newmapall,curmin]=convandthres(mapall,curmin,dict,CIflat,EBSDflat,K,active,dx,dy,dt,fid,rmspots,mexed);
    
k=1;
new=newmapall(:)~=mapall(:);
while k<=K
    mask1=(newmapall(:)==k);
    mask2=(mapall(:)==k);
    regsize=sum(mask1);
    change=sum(new(mask1|mask2));
    if change<3
        active(k)=0;
        k=k+1;
    else
        curmin(mask1)=fid*2;
        active(k)=1;
    changenum(k)=changenum(k)+change;
    num=changenum(k)/regsize;
    if regsize>10
    if (num)>=.4 && regsize>50
        changenum(k)=0;
        indices=find(mask1(:));
        if sum(CIflat(indices))>1e-4
            newind=datasample(indices,numsub,'Weights',CIflat(indices));
            EBSDtemp=EBSDflat(newind,:);
            if mexed
                [newg1, kap, ~] = VMFEMfast_mex(EBSDtemp, Pall,1,dict{k},kappa{k});
            else
                [newg1, kap, ~] = VMFEMfast(EBSDtemp, Pall,1,dict{k},kappa{k});
            end
            dict{k}=newg1;
            kappa{k}=kap;
        end
    end
    k=k+1;
    else
        if k<K
            dict{k}=dict{K};
            kappa{k}=kappa{K};
            mapall(mapall==k)=0;
            mapall(mapall==K)=k;
            newmapall(newmapall==k)=0;
            newmapall(newmapall==K)=k;
            changenum(k)=changenum(K);
            active(k)=active(K);
        end
        changenum(K)=[];
        kappa(K)=[];
        dict(K)=[];
        active(K)=[];
        K=K-1;
    end
    end
end
mapall=newmapall;
imagesc(mapall)
pause(1)
totalnum=sum(new(:));
if totalnum<2
    if dt<1/2^12
        break
    else
        dt=dt/2;
        active=ones(1,K);
    end
end

end