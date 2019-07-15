function [newmapall,newdict,newkappa,newK]=regionmerger(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,nt,numsub,betatol,conservative)
    if betatol>.75
        betatol=betatol/180*pi;
    end
    
if nargin<14
    conservative=.95;
end
    fac=1/(50*(dx+dy));
    K=size(dict,2);
    [neigharray,sizeneigh]=  findneigh(mapall,K);
    
    totneigh=sum(sizeneigh(:))/2;
    neighdist=zeros(totneigh,3);
    z=1;
    for k=1:K
        S=sizeneigh(k);
        for s=1:S
            cand=neigharray(k,s);
            if cand>k
               neighdist(z,:)=[k,cand,b2bmetric(dict(:,k),dict(:,cand))];
               z=z+1;
            end
        end
    end
    [~,I]=sort(neighdist(:,3));
    neighdist=neighdist(I,:);
    
    [M,N]=size(mapall);
    sqrtdt=sqrt(dt);
    w=ceil(fac*600*sqrtdt);
    [~,~,z]=size(EBSD);

    EBSDflat=E313toq(reshape(EBSD,[M*N,z])');
    
    
    [xbdcor,ybdcor,sizebdcor,coords,sizecoords,minmaxrowcol]=  bndcoords(mapall,K);%Gets coordinates
    numelemts=[K,ceil(8*fac*600*sqrt(dt)*sqrt(M*N/K))];


    fidK=zeros(numelemts);
    regK=zeros(numelemts);
    BindKsz=zeros(numelemts);
    SindKsz=zeros(numelemts);
    fidregKsz=zeros(K,1);
    bndsK=zeros(K,4);
    
    for k=1:K
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
    
    for k=1:K
        fullsz=fidregKsz(k);
        if fullsz>0 %solves the heat equation to do the convolution

            slinind=SindKsz(k,1:fullsz);
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
    map=1:K;
    active=ones(1,K);
    changed=zeros(1,K);
    for i=1:totneigh
        if neighdist(i,3)>betatol
            break
        end
        k1=neighdist(i,1);
        k2=neighdist(i,2);
        k1tok2=findbndenergy(k1,k2);
        k2tok1=findbndenergy(k2,k1);
        k1plusk2=findcombinedenergy(k1,k2);
        if k1plusk2<(k1tok2+k2tok1)*conservative
            root=findroot(map,k1);
            map(k2)=root;
            active(k2)=0;
            changed(root)=1;
        end
    end
    
    replacemap=zeros(1,K);
    newK=sum(active);
    newdict=zeros(4,newK);
    newkappa=zeros(1,newK);
    newchanged=zeros(1,newK);
    newk=1;
    for kz=1:K
        if active(kz)
            replacemap(kz)=newk;
            newdict(:,newk)=dict(:,kz);
            newkappa(newk)=kappa(kz);
            newchanged(newk)=changed(kz);
            newk=newk+1;
        else
            replacemap(kz)=replacemap(findroot(map,kz));
        end
    end
            
    newmapall=changemap(mapall,replacemap);
    [~,~,~,coords,sizecoords,~]=  bndcoords(newmapall,newK);
    for k=1:newK
        if newchanged(k)
            csize=sizecoords(k);
            indices=coords(k,1:csize);
            if sum(CI(indices))>1e-4
                [newg1,kap]=estimatebeta(csize,CI(indices),EBSDflat(:,indices),beta(indices),dict(:,k),kappa(k),numsub);
                newdict(:,k)=newg1;
                newkappa(k)=kap;
            end  
        end
    end
            
            
            
            
        
        
    function bndenergy=findbndenergy(k1,k2)
        fullsz=fidregKsz(k1);
        bndenergy=0;
        if fullsz>0
            slinind=SindKsz(k1,1:fullsz);
            miny=bndsK(k1,1);
            maxy=bndsK(k1,2);
            minx=bndsK(k1,3);
            maxx=bndsK(k1,4);
            m=maxy-miny+1;
            n=maxx-minx+1;
            indy=miny:maxy;
            indx=minx:maxx;
            endmask=(mapall(indy,indx)==k2);
            temp=endmask*0;
            temp(slinind)=1/sqrtdt*(regK(k1,1:fullsz));
            bndenergy=sum(temp(endmask));

            csize=sizecoords(k1);
            indices=coords(k1,1:csize);
            bmask=beta(indices);
            alphacoord=indices(~bmask);
            betacoord=indices(bmask);
            if ~isempty(alphacoord)
                bndenergy=bndenergy+sum(fid*CI(alphacoord).*alpbmetric(EBSDflat(:,alphacoord),dict(:,k1)));
            end
            if ~isempty(betacoord)
                bndenergy=bndenergy+sum(fid*CI(betacoord).*b2bmetric(EBSDflat(:,betacoord),dict(:,k1)));
            end
        end
    end 


    function combinedenergy=findcombinedenergy(k1,k2)
        csize1=sizecoords(k1);
        indices1=coords(k1,1:csize1);
        csize2=sizecoords(k2);
        indices2=coords(k2,1:csize2);
        csize=csize1+csize2;
        indices=[indices1 indices2];
        if sum(CI(indices))>1e-4
            [betacand,~]=estimatebeta(csize,CI(indices),EBSDflat(:,indices),beta(indices),dict(:,k1),kappa(k1),numsub);
        else
            betacand=dict(:,k1);
        end  
        combinedenergy=0;
        bmask=beta(indices);
        alphacoord=indices(~bmask);
        betacoord=indices(bmask);
        if ~isempty(alphacoord)
            combinedenergy=combinedenergy+sum(fid*CI(alphacoord).*alpbmetric(EBSDflat(:,alphacoord),betacand));
        end
        if ~isempty(betacoord)
            combinedenergy=combinedenergy+sum(fid*CI(betacoord).*b2bmetric(EBSDflat(:,betacoord),betacand));
        end

    end 

end
