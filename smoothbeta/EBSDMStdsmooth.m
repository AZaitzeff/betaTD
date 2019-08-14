function [mapall,smoothEBSD,energy,flag]=EBSDMStdsmooth(mapall,EBSD,CI,beta,smoothEBSD,K,fid,DT,dx,dy,dtstop,nt,between,ntsm,dtsm)
%This function is the main workhorse of the code
energy=inf;
flag=1;%if flag equals zero dt is too large
dts=2.^linspace(log2(DT),log2(dtstop),4);%limit at dt goes to 0



totalK=K;
current=ones(1,K,'logical');
[M,N]=size(mapall);
[~,~,z]=size(EBSD);
EBSDflat=E313toq(reshape(EBSD,[M*N,z])');

MAXITER=125;
fac=1/(50*(dx+dy));% variable 

[xbdcor,ybdcor,sizebdcor,coords,sizecoords,minmaxrowcol]=  bndcoords(mapall,K);%Gets coordinates
ratio=zeros(K,1);
for k=1:K
   boxsize=(minmaxrowcol(k,2)-minmaxrowcol(k,1)+1)*(minmaxrowcol(k,4)-minmaxrowcol(k,3)+1);
   ratio(k)=sizecoords(k)/boxsize; 
end
crdelemts=size(coords);
changecounter=zeros(K,1);
numelemts=[K,ceil(8*fac*600*sqrt(DT)*sqrt(M*N/K))];
curmin=zeros(M,N);


fidK=zeros(numelemts);
regK=zeros(numelemts);
BindKsz=zeros(numelemts);
SindKsz=zeros(numelemts);
fidregKsz=zeros(K,1);

bndsK=zeros(K,4);

for dt=dts
    
    sqrtdt=sqrt(dt);
    w=ceil(fac*600*sqrtdt);

    active=ones(1,K);
    activecounter=K;
    for k=1:K
        updatebnds(k);
    end
    for k=1:K
        updateFID(k,1);
    end
    
    for t=1:MAXITER
        %t
        curmin(:)=2*fid;
        for k=1:K
            if current(k)
                csize=sizecoords(k);
                indices=coords(k,1:csize);
                fullsz=fidregKsz(k);
                linind=BindKsz(k,1:fullsz);
                safe=setdiff(indices,linind);
                curmin(safe)=-inf;
            end

        end
        for times=1:between

            for k=1:K
                fullsz=fidregKsz(k);
                if active(k) && fullsz>0 %solves the heat equation to do the convolution

                    slinind=SindKsz(k,1:fullsz);
                    linind=BindKsz(k,1:fullsz);
                    miny=bndsK(k,1);
                    maxy=bndsK(k,2);
                    minx=bndsK(k,3);
                    maxx=bndsK(k,4);
                    m=maxy-miny+1;
                    n=maxx-minx+1;
                    indy=miny:maxy;
                    indx=minx:maxx;
                    smallu=1*(mapall(indy,indx)==k);
                    mask=smallu*0;
                    mask(slinind)=1;
                    [xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,m,n);
                    newu=TSz(smallu,dt,nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n,0);
                    regK(k,1:fullsz)=newu(slinind);
                end

            end

            for k=1:K %Thresholds
                if current(k)
                    fullsz=fidregKsz(k);
                    linind=BindKsz(k,1:fullsz);
                    for ind=1:fullsz
                        temp=2/sqrtdt*(-regK(k,ind))+fidK(k,ind);
                        canind=linind(ind);
                        if temp<curmin(canind)
                            curmin(canind)=temp;
                            mapall(canind)=k;
                        end
                    end
                    
                end
            end
            %imagesc(mapall);colorbar%for testing
            %pause(1)
        end



        [xbdcor,ybdcor,sizebdcor,coordsa,sizecoordsa,minmaxrowcol]=  bndcoords(mapall,K); 
        p=max(sizecoordsa);
        if p>crdelemts(2)
            crdelemts(2)=ceil(p*1.1);
            [coords,sizecoords]=growarray(coords,sizecoords,K,crdelemts,-1,1);
        end

        for k=1:K
            if current(k)
                csize=sizecoordsa(k);

                if csize>10


                    oldcsize=sizecoords(k);
                    oldindices=coords(k,1:oldcsize);
                    indices=coordsa(k,1:csize);
                    change=numel(setxor(indices,oldindices));
                    perchange=(change/csize);
                    if perchange<.005 || change<4 % see which regions have not changed (this is for speed alone)
                        activecounter=activecounter-active(k);
                        active(k)=0;

                    else

                        valcng=changecounter(k)+perchange;
                        changecounter(k)=valcng;
                        activecounter=activecounter+1-active(k);
                        active(k)=1;

                        if valcng>.2 %updates betas
                            changecounter(k)=0;
                            
                            if sum(CI(indices))>1e-4
                            [u]=estso3smoothonegrain(M,N,indices,csize,EBSDflat,beta,CI,smoothEBSD,ntsm,dtsm,fid);
                            smoothEBSD(:,indices)=u;
                                
                            
                            end  
                        end


                    end
                    sizecoords(k)=sizecoordsa(k);
                    coords(k,1:csize)=coordsa(k,1:csize);
                    boxsize=(minmaxrowcol(k,2)-minmaxrowcol(k,1)+1)*(minmaxrowcol(k,4)-minmaxrowcol(k,3)+1);
                    ratio(k)=sizecoords(k)/boxsize;
                else
                    current(k)=0;
                    activecounter=activecounter-active(k);
                    active(k)=0;
                    fidregKsz(k)=0;
                    sizebdcor(k)=0;
                    sizecoords(k)=0;
                    totalK=totalK-1;
                end
            end
        end
        if (K/totalK)>=1.5
            simplify();
        end
        if activecounter==0
            break
        end

        for k=1:K
            if active(k)
                updatebnds(k);
                updateFID(k,1); 
            end
        end
        
    if prctile(ratio(current),25)<.1 && t>1
        flag=0;
        return
    end
    
    end
    simplify();
    MAXITER=ceil(MAXITER/2);
end

% for k=1:K
%     csize=sizecoords(k);
%     indices=coords(k,1:csize);
%     if sum(CI(indices))>1e-4
%         [newg1,kap]=estimatebeta(csize,CI(indices),EBSDflat(:,indices),beta(indices),dict(:,k),kappa(k),numsub);
%         dict(:,k)=newg1;
%         kappa(k)=kap;
%     end  
% end

for k=1:K
    updatebnds(k);
    updateFID(k,0); 
end

% energy=0;
% for k=1:K
%     fullsz=fidregKsz(k);
%     if fullsz>0
%         slinind=SindKsz(k,1:fullsz);
%         miny=bndsK(k,1);
%         maxy=bndsK(k,2);
%         minx=bndsK(k,3);
%         maxx=bndsK(k,4);
%         m=maxy-miny+1;
%         n=maxx-minx+1;
%         indy=miny:maxy;
%         indx=minx:maxx;
%         endmask=(mapall(indy,indx)~=k);
%         temp=endmask*0;
%         temp(slinind)=1/sqrtdt*(regK(k,1:fullsz));
%         energy=energy+sum(temp(endmask));
%         csize=sizecoords(k);
%         indices=coords(k,1:csize);
%         bmask=beta(indices);
%         alphacoord=indices(~bmask);
%         betacoord=indices(bmask);
%         if ~isempty(alphacoord)
%             energy=energy+sum(fid*CI(alphacoord).*alpbmetric(EBSDflat(:,alphacoord),dict(:,k)));
%         end
%         if ~isempty(betacoord)
%             energy=energy+sum(fid*CI(betacoord).*b2bmetric(EBSDflat(:,betacoord),dict(:,k)));
%         end
%     end
% end



    function simplify()
        if K>totalK
        map=zeros(1,K);
        newk=1;
        
        for kz=1:K
            if current(kz)
                
                changecounter(newk)=changecounter(kz);
                active(newk)=active(kz);
                
                sc=sizecoords(kz);
                sizecoords(newk)=sc;
                coords(newk,1:sc)=coords(kz,1:sc);
                
                num=fidregKsz(kz);
                fidregKsz(newk)=num;
                fidK(newk,1:num)=fidK(kz,1:num);
                regK(newk,1:num)=regK(kz,1:num);
                BindKsz(newk,1:num)=BindKsz(kz,1:num);
                SindKsz(newk,1:num)=SindKsz(kz,1:num);
                
                sbd=sizebdcor(kz);
                sizebdcor(newk)=sbd;
                xbdcor(newk,1:sbd)=xbdcor(kz,1:sbd);
                ybdcor(newk,1:sbd)=ybdcor(kz,1:sbd);
                minmaxrowcol(newk,:)=minmaxrowcol(kz,:);
                
                bndsK(newk,:)=bndsK(kz,:);
                map(kz)=newk;
                ratio(newk)=ratio(kz);
                
                
                newk=newk+1;
                
            end
        end
        mapall=changemap(mapall,map);
        changecounter(newk:K)=[];
        active(newk:K)=[];

        sizecoords(newk:K)=[];
        coords(newk:K,:)=[];

        fidregKsz(newk:K)=[];
        fidK(newk:K,:)=[];
        regK(newk:K,:)=[];
        BindKsz(newk:K,:)=[];
        SindKsz(newk:K,:)=[];

        sizebdcor(newk:K)=[];
        xbdcor(newk:K,:)=[];
        ybdcor(newk:K,:)=[];
        minmaxrowcol(newk:K,:)=[];

        bndsK(newk:K,:)=[];
        ratio(newk:K)=[];
        
        
        K=totalK;
        current=ones(1,K,'logical');
        end
    end
    function updateFID(k,smooth)
        fullsz=fidregKsz(k);
        if fullsz>0
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
            bmask=beta(linind);
            alphacoord=linind(~bmask);
            betacoord=linind(bmask);
            if ~isempty(alphacoord)
                slin=slinind(~bmask);
                nfid=sum(~bmask);
                for ifid=1:nfid
                    mask(slin(ifid))=fid*CI(alphacoord(ifid)).*alpbmetric(EBSDflat(:,alphacoord(ifid)),smoothEBSD(:,alphacoord(ifid)));
                end
            end
            if ~isempty(betacoord)
                slin=slinind(bmask);
                nfid=sum(bmask);
                for ifid=1:nfid
                    mask(slin(ifid))=fid*CI(betacoord(ifid)).*alpbmetric(EBSDflat(:,betacoord(ifid)),smoothEBSD(:,betacoord(ifid)));
                end
            end
            if smooth
                newu=TSz(mask,min(dt/16,2^-12),nt,dx,dy,xdir,ydir,xsizes,ysizes,m,n,1);
            else
                newu=mask;
            end
            fidK(k,1:fullsz)=newu(slinind);
        end
    end
    function updatebnds(k)
        total=sizebdcor(k);
        if total>0
            [linind,slinind,fullsz,bnds]=findboundary(w,minmaxrowcol(k,:),xbdcor(k,1:total)',ybdcor(k,1:total)',M,N);
            if fullsz>numelemts(2)
                numelemts(2)=ceil(fullsz*1.1);
                [BindKsz,fidregKsz]=growarray(BindKsz,fidregKsz,K,numelemts,-1,1);
                [SindKsz,fidregKsz]=growarray(SindKsz,fidregKsz,K,numelemts,-1,1);
                [fidK,fidregKsz]=growarray(fidK,fidregKsz,K,numelemts,-1,1);
                [regK,fidregKsz]=growarray(regK,fidregKsz,K,numelemts,-1,1);
            end
            [sortedind,order]=sort(linind);
            BindKsz(k,1:fullsz)=sortedind;
            SindKsz(k,1:fullsz)=slinind(order);
            fidregKsz(k)=fullsz;
            bndsK(k,:)=bnds;
        else
            bndsK(k,:)=[0,0,0,0];
            fidregKsz(k)=0;
        end
    end
end