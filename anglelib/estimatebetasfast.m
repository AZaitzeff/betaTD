function [dict,kappa]=estimatebetasfast(EBSD,CI,map,max_init,mexed)
    numsub=200;
    if nargin<4
        max_init=8;
    end
    if nargin<5
        mexed=0;
    end
    K = max(map(:));
    dict=cell(1,K);
    kappa=cell(1,K);

[m,n,z]=size(EBSD);
EBSDflat=reshape(EBSD, [m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI, [m*n,1]);


T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
for z=1:K
    indices=find(z==map);
    if numel(indices)>10
        CItemp=CIflat(indices);
        if sum(CIflat(indices))==0
            CItemp=CItemp+1;
        end
        newind=datasample(indices,numsub,'Weights',CItemp);
        EBSDtemp=EBSDflat(newind,:);
        if mexed
            [mu, kap, ~] = VMFEMfast_mex(EBSDtemp, Pall,max_init,[0,0,0,0],1);
        else
            [mu, kap, ~] = VMFEMfast(EBSDtemp, Pall,max_init,[0,0,0,0],1);
        end
        dict{z}=mu;
        kappa{z}=kap;
    else
        dict{z}=[1,0,0,0];
        kappa{z}=1;
    end
end