function [dict,kappa,IQmean,betamean]=estimatemeanfast(EBSDflat,CI,IQ,beta,map,K,numsub)

    max_init=4;

    dict=zeros(4,K);
    kappa=zeros(1,K);
    IQmean=zeros(1,K);
    betamean=zeros(1,K,'logical');


Pm1=getsymmetries('hex');
Pm2=getsymmetries('cubic');

[~,~,~,coords,sizecoords,~]=  bndcoords(map,K);
for k=1:K
    indices=coords(k,1:sizecoords(k));
    if numel(indices)>10
        CItemp=CI(indices);
        if sum(CI(indices))==0
            CItemp=CItemp+1;
        end
        newind=datasamplez(indices,numsub,CItemp);
        EBSDtemp=EBSDflat(:,newind);
        mask=beta(newind);
        alphaEBSD=EBSDtemp(:,~mask);
        betaEBSD=EBSDtemp(:,mask);
        [mu, kap, ~] = VMFEMzfast(alphaEBSD, Pm1,betaEBSD, Pm2,max_init,[1;0;0;0],1);
        %misoriena=alpbmetric(alphaEBSD,mu);
        %misorienb=b2bmetric(betaEBSD,mu);
        %tol=1.5*median([misoriena misorienb]);
        %maska=misoriena<tol;
        %maskb=misorienb<tol;
        %[mu, kap, ~] = VMFEMzfast(alphaEBSD(:,maska), Pall,betaEBSD(:,maskb), Pm,1,mu,kap);
        dict(:,k)=mu;
        kappa(k)=kap;
        IQmean(k)=mean(IQ(newind));
        betamean(k)=mean(mask)>.5;
        
    else
        dict(:,k)=[1,0,0,0];
        kappa(k)=50;
        IQmean(k)=0;
        betamean(k)=0;
    end
end
end
