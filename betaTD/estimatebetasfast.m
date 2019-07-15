function [dict,kappa]=estimatebetasfast(EBSDflat,CIflat,beta,map,K,numsub)
    
   %{
    Given regions finds the beta mean of each one
   %}

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
    N=sizecoords(k);
    indices=coords(k,1:N);
    
    if N>10
        CItemp=CIflat(indices);
        if sum(CIflat(indices))==0
            CItemp=CItemp+1;
        end
        [mu,kap]=estimatebeta(N,CItemp,EBSDflat(:,indices),beta(indices),[],[],numsub);
        
        dict(:,k)=mu;
        kappa(k)=kap;
    else
        dict(:,k)=[1,0,0,0];
        kappa(k)=1;
    end
end
