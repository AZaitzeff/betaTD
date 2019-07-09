function [mapall,dict,kappa]=EBSDMStdz(mapall,EBSD,CI,dict,kappa,fid,dx,dy,DT,rmspots,dtstop,mexed,numsub)
if nargin<8
    dx=1/100;
    dy=1/100;
end
if nargin<9
    DT=.02;
end
if nargin<13
    numsub=200;
end
if nargin<12
    mexed=0;
end
if nargin<11
    dtstop=2^(-12);
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
[M,N]=size(mapall);%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[M*N,1]);
MAXITER=1000;
%curmin=ones(M,N)*fid*2+1;
active=ones(1,K);
[xbdcor,ybdcor,sizebdcor,coordsb,sizecorb,~,minmaxrowcol]=  bndcoords(mapall,K);
coordsstill=coordsb;
sizecorstill=sizecorb;
fidterm=ones(K,M*N)*inf;

ene=zeros(K,K);
for k1=1:K
    for k2=k1+1:K
        energy=readshockley(dict{k1},dict{k2});
        ene(k1,k2)=energy;
        ene(k2,k1)=energy;
    end
end
for t=1:MAXITER
    

fac=1/(50*(dx+dy));
w=ceil(fac*400*sqrt(dt));
oldmapall=mapall;
vals=ones(K,M*N)*inf;
for k=1:K
    if active(k)
        total=sizebdcor(k);

        [xdir,ydir,smallu,linind,slinind,m,n,xsizes,ysizes]=findboundary(oldmapall,k,w,minmaxrowcol(k,:),xbdcor(k,1:total)',ybdcor(k,1:total)',M,N,rmspots,mexed);
        if m>0
            newu=ADIz(smallu*1,dt,dx,dy,xdir,ydir,xsizes,ysizes,m,n);
            %mask=(newu(slinind)>.025 & newu(slinind)<.995);
            %slinind=slinind(mask);
            %linind=linind(mask);
            mask=isinf(fidterm(k,linind));
            fidterm(k,linind(mask))=CIflat(linind(mask)).*alpbmetric(EBSDflat(linind(mask),:),dict{k})';
            %fidterm{k}=temp;
            mask=isinf(vals(k,linind));
            vals(k,linind(mask))=0;
            vals(k,linind)=vals(k,linind)+fid*fidterm(k,linind);
            for kn=1:K
                mask=isinf(vals(kn,linind));
                vals(kn,linind(mask))=0;
                vals(kn,linind)=vals(kn,linind)+ene(k,kn)*2/sqrt(dt)*(newu(slinind))';
            end
            %[minval,I]=findminz(curmin(linind),S);
            %curmin(linind(I))=minval(I);
            %mapall(linind(I))=k;
        end
    end
end
[~,I]=min(vals,[],1);
mapall=reshape(I,[M N]);
imagesc(mapall)
pause(1)
[xbdcor,ybdcor,sizebdcor,coordsa,sizecora,neighbors,minmaxrowcol]=  bndcoords(mapall,K);

k=1;
totalnum=0;
while k<=K
    regsize=sizecora(k);
    indices=coordsa(k,1:regsize);
    change=numel(setxor(indices,coordsb(k,1:sizecorb(k))));
    totalnum=totalnum+change;
    if change<3 && regsize>10
        %active(k)=0;
        k=k+1;
    else
        %curmin(indices)=fid*2+1;
        active(k)=1;
    
    if regsize>10
        num=numel(setxor(indices,coordsstill(k,1:sizecorstill(k))))/regsize;
    if (num)>=.4 && regsize>50
        coordsstill(k,:)=coordsa(k,:);
        sizecorstill(k)=sizecora(k);
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
            fidterm(k,:)=inf;
            for k2=1:K
                if k~=k2
                    energy=readshockley(dict{k},dict{k2});
                    ene(k,k2)=energy;
                    ene(k2,k)=energy;
                end

            end

        end
    end
    k=k+1;
    else
        if k<K
            dict{k}=dict{K};
            kappa{k}=kappa{K};
            mapall(mapall==k)=0;
            mapall(mapall==K)=k;
            active(k)=active(K);
            ene(k,:)=ene(K,:);
            ene(:,k)=ene(:,k);
            fidterm(k,:)=fidterm(K,:);
            coordsstill(k,:)=coordsstill(K,:);
            coordsa(k,:)=coordsa(K,:);
            coordsb(k,:)=coordsb(K,:);
            sizecorstill(k)=sizecorstill(K);
            sizecora(k)=sizecora(K);
            sizecorb(k)=sizecorb(K);
        end
        kappa(K)=[];
        dict(K)=[];
        active(K)=[];
        fidterm(K,:)=[];
        coordsstill(K,:)=[];
        coordsa(K,:)=[];
        coordsb(K,:)=[];
        sizecorstill(K)=[];
        sizecora(K)=[];
        sizecorb(K)=[];
        K=K-1;
    end
    end
end
coordsb=coordsa;
sizecorb=sizecora;

if totalnum<2
    if dt<dtstop
        break
    else
        dt=dt/2;
        active=ones(1,K);
    end
end

end