function [dict,kappa]=estimatebetasfast(EBSDflat,CIflat,beta,map,K,numsub)

    max_init=8;

    dict=zeros(4,K);
    kappa=zeros(1,K);



T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=Pm(:,:,ceil(i/6))*T(:,:,mod(i-1,6)+1)';
end
[~,~,~,coords,sizecoords,~]=  bndcoords(map,K);
for k=1:K
    indices=coords(k,1:sizecoords(k));
    if numel(indices)>10
        CItemp=CIflat(indices);
        if sum(CIflat(indices))==0
            CItemp=CItemp+1;
        end
        newind=datasamplez(indices,numsub,CItemp);
        EBSDtemp=EBSDflat(:,newind);
        mask=beta(newind);
        alphaEBSD=EBSDtemp(:,~mask);
        betaEBSD=EBSDtemp(:,mask);
        [mu, kap, ~] = VMFEMzfast(alphaEBSD, Pall,betaEBSD, Pm,max_init,[1;0;0;0],1);
        %misoriena=alpbmetric(alphaEBSD,mu);
        %misorienb=b2bmetric(betaEBSD,mu);
        %tol=1.5*median([misoriena misorienb]);
        %maska=misoriena<tol;
        %maskb=misorienb<tol;
        %[mu, kap, ~] = VMFEMzfast(alphaEBSD(:,maska), Pall,betaEBSD(:,maskb), Pm,1,mu,kap);
        dict(:,k)=mu;
        kappa(k)=kap;
    else
        dict(:,k)=[1,0,0,0];
        kappa(k)=1;
    end
end
