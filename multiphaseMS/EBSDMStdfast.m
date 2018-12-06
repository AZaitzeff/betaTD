function [mapall,newdict,newkappa,energy]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,DT,dx,dy,dtstop,nt,numsub)


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
    dtstop=2^(-10);
end
between=10;
dt=DT;
%mult=ceil(log2(DT/dtstop));
%fid=fid/2^(mult/2);
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
energy=0;
K=max(mapall(:));
newmapall=mapall;
current=ones(1,K);
[M,N]=size(mapall);
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,z]);
EBSDflat=E313toq(EBSDflat);
MAXITER=200;
%curmin=ones(M,N)*fid*2;
active=ones(1,K);
activecounter=K;
fac=1/(50*(dx+dy));
w=ceil(fac*600*sqrt(dt));
[xbdcor,ybdcor,sizebdcor,coords,sizecoords,minmaxrowcol]=  bndcoords(mapall,K);
%bdelemts=size(ybdcor);
crdelemts=size(coords);
changecounter=zeros(K,1);
%numelemts=[K,ceil(M*N)];
numelemts=[K,ceil(8*w*sqrt(M*N/K))];
curmin=zeros(M,N);
fidK=zeros(numelemts);
regK=zeros(numelemts);
BindKsz=zeros(numelemts);
SindKsz=zeros(numelemts);
fidregKsz=zeros(K,1);
bndsK=zeros(K,4);

for k=1:K
    updatebnds(k);
end
for k=1:K
    updateFID(k);
end
for t=1:MAXITER
    
for times=1:between
    newmapall(:)=mapall(:);
    curmin(:)=fid*2;
    for k=1:K
        if current(k)
            fullsz=fidregKsz(k);
            slinind=SindKsz(k,1:fullsz);
            linind=BindKsz(k,1:fullsz);
            if active(k)
                miny=bndsK(k,1);
                maxy=bndsK(k,2);
                minx=bndsK(k,3);
                maxx=bndsK(k,4);
                m=maxy-miny+1;
                n=maxx-minx+1;
                indy=miny:maxy;
                indx=minx:maxx;
                smallu=1*(newmapall(indy,indx)==k);
                mask=smallu*0;
                mask(slinind)=1;
                [xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,m,n);
                newu=TSz(smallu,dt,nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n,0);
                regK(k,1:fullsz)=newu(slinind);
            end

            S=2/sqrt(dt)*(-regK(k,1:fullsz))+fidK(k,1:fullsz);
            [minval,I]=findminz(curmin(linind),S);
            curmin(linind(I))=minval(I);
            mapall(linind(I))=k;
        end
    end
    %imagesc(mapall)
    %pause(1)
end

[xbdcor,ybdcor,sizebdcor,coordsa,sizecoordsa,minmaxrowcol]=  bndcoords(mapall,K); 
p=max(sizecoordsa);
if p>crdelemts(2)
    crdelemts(2)=ceil(p*1.1);
    [coords,sizecoords]=growarray(coords,sizecoords,K,crdelemts,-1,1);
end
% p=max(sizebdcora);
% if p>bdelemts(2)
%     bdelemts(2)=ceil(p*1.1);
%     [ybdcora,sizebdcor]=growarray(ybdcora,sizebdcor,K,bdelemts,-1,1);
%     [xbdcora,sizebdcor]=growarray(xbdcora,sizebdcor,K,bdelemts,-1,1);
% end
for k=1:K
    if current(k)
        
        csize=sizecoordsa(k);

        if csize>10
            
            
            oldcsize=sizecoords(k);
            oldindices=coords(k,1:oldcsize);
            indices=coordsa(k,1:csize);
            change=numel(setxor(indices,oldindices));
            perchange=(change/csize);
            if perchange<.01 || change<10
                activecounter=activecounter-active(k);
                active(k)=0;
                
            else
                
                valcng=changecounter(k)+perchange;
                changecounter(k)=valcng;
                activecounter=activecounter+1-active(k);
                active(k)=1;
                    
                if valcng>.2
                    changecounter(k)=0;
                    if sum(CI(indices))>1e-4
                        newind=datasamplez(indices,numsub,CI(indices));
                        EBSDtemp=EBSDflat(newind,:);
                        mask=beta(newind);
                        alphaEBSD=EBSDtemp(~mask,:);
                        betaEBSD=EBSDtemp(mask,:);
                        [newg1, kap, ~] = VMFEMzfast(alphaEBSD, Pall,betaEBSD, Pm,1,dict(k,:),kappa(k));
                        %[newg1, kap, ~] = VMFEMfast(EBSDtemp, Pall,1,dict(k,:),kappa(k));
                        dict(k,:)=newg1;
                        kappa(k)=kap;
                    end  
                end
                       

            end
            sizecoords(k)=sizecoordsa(k);
            coords(k,1:csize)=coordsa(k,1:csize);
        else
            current(k)=0;
            activecounter=activecounter-active(k);
            active(k)=0;
            fidregKsz(k)=0;
            sizebdcor(k)=0;
            sizecoords(k)=0;
        end
    end
end
for k=1:K
    if active(k)
        updatebnds(k);
        updateFID(k); 
    end
end

if activecounter==0
    if dt<=dtstop
        energy=0;
        for k=1:K
            if current(k)
                fullsz=fidregKsz(k);
                slinind=SindKsz(k,1:fullsz);
                miny=bndsK(k,1);
                maxy=bndsK(k,2);
                minx=bndsK(k,3);
                maxx=bndsK(k,4);
                m=maxy-miny+1;
                n=maxx-minx+1;
                indy=miny:maxy;
                indx=minx:maxx;
                endmask=(newmapall(indy,indx)~=k);
                temp=endmask*0;
                temp(slinind)=1/sqrt(dt)*(regK(k,1:fullsz));
                energy=energy+sum(temp(endmask));
                csize=sizecoordsa(k);
                indices=coordsa(k,1:csize);
                mask=beta(linind);
                energy=energy+sum(fid*CI(indices(~mask)).*alpbmetric(EBSDflat(indices(~mask),:),dict(k,:)));
                energy=energy+sum(fid*CI(indices(mask)).*b2bmetric(EBSDflat(indices(mask),:),dict(k,:)));
            end
        end
        break
    else
        %fid=fid*sqrt(2);
        dt=dt/2;
        w=ceil(fac*600*sqrt(dt));
        
        
        active=current;
        activecounter=sum(active);
        for k=1:K
            if current(k)
                updatebnds(k);
            end
        end
        for k=1:K
            if current(k)
                updateFID(k);
            end
        end
    end
end

end

if t==MAXITER
   energy=inf; 
end
newK=sum(current);
newdict=zeros(newK,4);
newkappa=zeros(newK,1);
map=1:K;
newk=1;
for k=1:K
    if current(k)
        map(k)=newk;
        newdict(newk,:)=dict(k,:);
        newkappa(newk)=kappa(k);
        newk=newk+1;
    end
end
mapall=changemap(mapall,map);

    function updateFID(k)
        fullsz=fidregKsz(k);
        slinind=SindKsz(k,1:fullsz);
        linind=BindKsz(k,1:fullsz);
        miny=bndsK(k,1);
        maxy=bndsK(k,2);
        minx=bndsK(k,3);
        maxx=bndsK(k,4);
        m=maxy-miny+1;
        n=maxx-minx+1;
        mask=zeros(m,n);
        mask(slinind)=1;
        [xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,m,n);
        mask=beta(linind);
        mask(slinind(~mask))=fid*CI(linind(~mask)).*alpbmetric(EBSDflat(linind(~mask),:),dict(k,:));
        mask(slinind(mask))=fid*CI(linind(mask)).*b2bmetric(EBSDflat(linind(mask),:),dict(k,:));
        newu=TSz(mask,dt/16,nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n,1);
        fidK(k,1:fullsz)=newu(slinind);
    end
    function updatebnds(k)
        total=sizebdcor(k);
        [linind,slinind,fullsz,bnds]=findboundary(w,minmaxrowcol(k,:),xbdcor(k,1:total)',ybdcor(k,1:total)',M,N);
        if fullsz>numelemts(2)
            numelemts(2)=ceil(fullsz*1.1);
            [BindKsz,fidregKsz]=growarray(BindKsz,fidregKsz,K,numelemts,-1,1);
            [SindKsz,fidregKsz]=growarray(SindKsz,fidregKsz,K,numelemts,-1,1);
            [fidK,fidregKsz]=growarray(fidK,fidregKsz,K,numelemts,-1,1);
            [regK,fidregKsz]=growarray(regK,fidregKsz,K,numelemts,-1,1);
        end
        BindKsz(k,1:fullsz)=linind;
        SindKsz(k,1:fullsz)=slinind;
        fidregKsz(k)=fullsz;
        bndsK(k,:)=bnds;
        
    end
end