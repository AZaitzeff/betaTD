function [mapall,dict,kappa]=EBSDMStdfast(mapall,EBSD,CI,dict,kappa,fid,DT,dx,dy,dtstop,nt,numsub)

%TO DO: store fidelity cals in a smaller matrix. map to lower indices. see
%if boundary stays within a certian region

if nargin<8
    dx=1/100;
    dy=1/100;
end
if nargin<7
    DT=.02;
end
if nargin<11
    nt=5;
end
if nargin<12
    numsub=400;
end
if nargin<10
    dtstop=2^(-12);
end
dt=DT;
%mult=ceil(log2(DT/dtstop));
%fid=fid/2^(mult/2);
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
MAXITER=1000;
curmin=ones(m,n)*fid*2;
active=ones(1,K);
[xbdcor,ybdcor,sizebdcor,coords1,sizecor1,minmaxrowcol]=  bndcoords(mapall,K);
coordsfix=coords1;
sizecorfix=sizecor1;
[M,N]=size(mapall);
fac=1/(50*(dx+dy));
w=ceil(fac*400*sqrt(dt));
%w1=ceil(fac*400*sqrt(dt));
%w2=ceil(fac*800*sqrt(dt));

for t=1:MAXITER

newmapall=mapall;
for k=1:K
    
    if active(k)
        total=sizebdcor(k);
        if total>5
            [xdir,ydir,xsizes,ysizes,smallu,linind,slinind,m,n]=findboundary(newmapall,k,w,minmaxrowcol(k,:),xbdcor(k,1:total)',ybdcor(k,1:total)',M,N);
        
        
        
            newu=TSz(smallu*1,dt,nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n);
            %mask=(newu(slinind)>.025 & newu(slinind)<.995);
            %slinind=slinind(mask);
            %linind=linind(mask);
            S=2/sqrt(dt)*(-newu(slinind))+fid*CI(linind).*alpbmetric(EBSDflat(linind,:),dict(k,:));
            [minval,I]=findminz(curmin(linind),S);
            curmin(linind(I))=minval(I);
            mapall(linind(I))=k;
        end
    end
end
% imagesc(mapall)
% pause(1)
[xbdcor,ybdcor,sizebdcor,coords2,sizecor2,minmaxrowcol]=  bndcoords(mapall,K);    
k=1;
totalnum=0;
while k<=K
    regsize=sizecor2(k);
    indices=coords2(k,1:regsize);
    change=numel(setxor(indices,coords1(k,1:sizecor1(k))));
    totalnum=totalnum+change;
    
    if change<3 && regsize>10
        active(k)=0;
        k=k+1;
    else
        curmin(indices)=fid*2;
        active(k)=1;

    if regsize>10
        num=numel(setxor(indices,coordsfix(k,1:sizecorfix(k))))/regsize;
    if (num)>=.2 && regsize>50
        if sum(CI(indices))>1e-4
            newind=datasamplez(indices,numsub,CI(indices));
            EBSDtemp=EBSDflat(newind,:);
            [newg1, kap, ~] = VMFEMfast(EBSDtemp, Pall,1,dict(k,:),kappa(k));
            dict(k,:)=newg1;
            kappa(k)=kap;
        end
    end
    k=k+1;
    else
        if k<K
            dict(k,:)=dict(K,:);
            kappa(k)=kappa(K);
            mapall(mapall==k)=0;
            mapall(mapall==K)=k;
            active(k)=active(K);
            coordsfix(k,:)=coordsfix(K,:);
            coords1(k,:)=coords1(K,:);
            coords2(k,:)=coords2(K,:);
            sizecorfix(k)=sizecorfix(K);
            sizecor1(k)=sizecor1(K);
            sizecor2(k)=sizecor2(K);
            xbdcor(k,:)=xbdcor(K,:);
            ybdcor(k,:)=ybdcor(K,:);
            sizebdcor(k,:)=sizebdcor(K,:);
            minmaxrowcol(k,:)=minmaxrowcol(K,:);
        end
        kappa(K)=[];
        dict(K,:)=[];
        active(K)=[];
        coordsfix(K,:)=[];
        coords1(K,:)=[];
        coords2(K,:)=[];
        sizecorfix(K)=[];
        sizecor1(K)=[];
        sizecor2(K)=[];
        xbdcor(K,:)=[];
        ybdcor(K,:)=[];
        sizebdcor(K,:)=[];
        minmaxrowcol(K,:)=[];
        K=K-1;
    end
    end
end
coords1=coords2;
sizecor1=sizecor2;


if totalnum<2
    if dt<dtstop
        break
    else
        %fid=fid*sqrt(2);
        dt=dt/2;
        w=ceil(fac*600*sqrt(dt));
        active=ones(1,K);
    end
end

end