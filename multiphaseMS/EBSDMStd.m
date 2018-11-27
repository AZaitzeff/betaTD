function [mapall,dict,kappa]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,DT,dx,dy,dtstop,enec,numsub)
%vidObj = VideoWriter(['grim'], 'MPEG-4');
%vidObj.Quality = 100;
%vidObj.FrameRate = 4;
%open(vidObj);
%imagesc(mapall)
%writeVideo(vidObj, getframe(gca));
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
    enec=0;
end
if nargin<12
    numsub=400;
end
flag=1;
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
    S{k}(:)=CIflat(:).*alpbmetric(EBSDflat(:,:),dict(k,:))';
end
S=td2dz(S,dt/4,dx,dy);

changeu=u;
lastu=u;
for t=1:MAXITER


%for k=1:K
%    mask=(ls{k}<.995)&(ls{k}>.01)&((S{k}>.95)|(S{k}<-.01));
%    S{k}(mask)=CIflat(mask(:)).*alpbmetric(EBSDflat(mask(:),:),dict(k,:))';
%end
for k=1:K
    mapall(u{k}>0)=k;
end

if enec
     ene=zeros(K,K);
     for k1=1:K
         for k2=(k1+1):K
             energy=readshockley(dict(k1,:),dict(k2,:));
             ene(k1,k2)=energy;
             ene(k2,k1)=energy;
         end
     end
end

[ls]=td2dz(u,dt,dx,dy);
phi=zeros(m,n,K);

for k=1:K
    phi(:,:,k)=fid*S{k};
    if ~enec
    phi(:,:,k)=phi(:,:,k)-2/(sqrt(dt))*ls{k};
    else
        
     for k2=1:K
         phi(:,:,k)=phi(:,:,k)+2/(sqrt(dt))*ene(k,k2)*ls{k2};
     end
    end
end

[~,argmin]=min(phi,[],3);
for k=1:K
    u{k}=argmin==k;
    
end

mapall=reshape(argmin, [m,n]);
%energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dt,dx,dy)
%imagesc(mapall)
%pause(1)
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
            newind=datasamplez(indices,numsub,CI(indices));
            EBSDtemp=EBSDflat(newind,:);
            [newg1, kap, ~] = VMFEMfast(EBSDtemp, Pall,1,dict(k,:),kappa(k));
            dict(k,:)=newg1;
            kappa(k)=kap;
            S{k}(:)=CIflat(:).*alpbmetric(EBSDflat(:,:),dict(k,:))';
            ls=td2dz({S{k}},dt/4,dx,dy);
            S{k}=ls{1};
        end
    end
    k=k+1;
    else
        dict(k,:)=dict(K,:);
        u{k}=u{K};
        kappa(k)=kappa(K);
        kappa(K)=[];
        u(K)=[];
        dict(K,:)=[];
        changeu{k}=changeu{K};
        lastu{k}=lastu{K};
        S{k}=S{K};
        S(K)=[];
        K=K-1;
    end
end

%title(['energy: ' num2str(energy) ' dt: ' num2str(dt)])
%writeVideo(vidObj, getframe(gca));
if totalnum<2
    if dt<=dtstop
        break
    elseif dt<=(dtstop*2)
        dt=dt/2;
        for k=1:K
            mapall(u{k}>0)=k;
        end
        [mapall,dict,kappa,K1]=paircheck(mapall,dict,kappa,CIflat,EBSDflat,K,dt,dx,dy,enec);
        K=K1;
        u=cell(1,K);
        S=cell(1,K);
        for k=1:K
            S{k}(:)=CIflat(:).*alpbmetric(EBSDflat(:,:),dict(k,:))';
            u{k}=mapall==k;
        end
    else
        dt=dt/2;
        for k=1:K
            S{k}(:)=CIflat(:).*alpbmetric(EBSDflat(:,:),dict(k,:))';
        end
        S=td2dz(S,dt/4,dx,dy);
    end
end
lastu=u;


end

for k=1:K
    mapall(u{k}>0)=k;
end
%close(vidObj);