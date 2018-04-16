function [dict]=estimatebetas(EBSD,CI,betas,map,update,dict,sub,numsub)
if numel(update)==0
    dict = containers.Map('KeyType','int32','ValueType','any');
    [clusterlist,~,~] = unique(map);
else
    clusterlist=update';
end
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
    CItemp=CIflat(indices);
    if sum(CIflat(indices))==0
        CItemp=CItemp+1;
    end
    if sub
        
        newind=datasample(indices,numsub,'Weights',CItemp);
        CItemp=CItemp*0+1;
    else
        newind=indices;
    end
    EBSDtemp=EBSDflat(newind,:);
    maskalpha=alphamask(newind);
    maskbeta=betamask(newind);
    if sum(maskbeta)==0
        [mu, ~, ~, ~] = VMFEM(EBSDtemp, Pall,CItemp,1,16);
    elseif sum(maskalpha)==0
        [mu, ~, ~, ~] = VMFEM(EBSDtemp, Pm,CItemp,1,16);
    else
        [mu, ~, ~, ~] = VMFEMz(EBSDtemp(maskalpha,:), Pall,CItemp(maskalpha),...
            EBSDtemp(maskbeta,:), Pm,CItemp(maskbeta),1,16);
    end
    dict(z)=mu;
end