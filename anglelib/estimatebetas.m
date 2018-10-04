function [dict,kappa]=estimatebetas(EBSD,CI,betas,map,update,dict,sub,numsub,max_init)
    
    if nargin<5
    update=[];
    dict=0;
    sub=1;
    end
    if nargin<8
        numsub=200;
    end
    if nargin<9
        max_init=4;
    end
if numel(update)==0
    dict = containers.Map('KeyType','int32','ValueType','any');

    [clusterlist,~,~] = unique(map);
else
    clusterlist=update';
end
kappa = containers.Map('KeyType','int32','ValueType','any');
[m,n,z]=size(EBSD);
EBSDflat=reshape(EBSD, [m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI, [m*n,1]);
betasflat=reshape(betas,[m*n,1]);
betamask=betasflat==1;
alphamask=betasflat==0;


T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
for z=clusterlist'
    indices=find(z==map);
    if numel(indices)>10
        CItemp=CIflat(indices);
        if sum(CIflat(indices))==0
            CItemp=CItemp+1;
        end
        if sub && numel(indices)>numsub
            newind=datasample(indices,numsub,'Weights',CItemp);
            CItemp=ones(numsub,1);
        else
            newind=indices;
        end
        EBSDtemp=EBSDflat(newind,:);
        maskalpha=alphamask(newind);
        maskbeta=betamask(newind);
        if sum(maskbeta)==0
            [mu, kap, ~, ~] = VMFEM(EBSDtemp, Pall,CItemp,1,max_init);
        elseif sum(maskalpha)==0
            [mu, kap, ~, ~] = VMFEM(EBSDtemp, Pm,CItemp,1,max_init);
        else
            [mu, kap, ~, ~] = VMFEMz(EBSDtemp(maskalpha,:), Pall,CItemp(maskalpha),...
                EBSDtemp(maskbeta,:), Pm,CItemp(maskbeta),1,max_init);
        end
        dict(z)=mu;
        kappa(z)=kap;
    else
        dict(z)=[1,0,0,0];
        kappa(z)=1;
    end
end