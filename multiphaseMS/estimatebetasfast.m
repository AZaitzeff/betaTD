function [dict,kappa]=estimatebetasfast(EBSDflat,CIflat,beta,map,K)

    max_init=16;

    numsub=400;

    dict=zeros(K,4);
    kappa=zeros(K,1);



T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
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
        EBSDtemp=EBSDflat(newind,:);
        mask=beta(newind);
        alphaEBSD=EBSDtemp(~mask,:);
        betaEBSD=EBSDtemp(mask,:);
        [mu, kap, ~] = VMFEMzfast(alphaEBSD, Pall,betaEBSD, Pm,max_init,[0,0,0,0],1);
        dict(k,:)=mu;
        kappa(k)=kap;
    else
        dict(k,:)=[1,0,0,0];
        kappa(k)=1;
    end
end
