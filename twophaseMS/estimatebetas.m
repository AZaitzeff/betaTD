function [dict]=estimatebetas(EBSD,CI,map,update,dict,numsub)
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


T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
for z=clusterlist'
    indices=find(z==map);
    if numsub<numel(indices)*3/4
        newind=datasample(indices,numsub,'Weights',CIflat(indices));
        CItemp=ones(size(CIflat(newind)));
    else
        newind=indices;
        CItemp=CIflat(newind);
    end
    [mu, ~, ~, ~] = VMFEM(EBSDflat(newind,:), Pall,CItemp,1,16);
    dict(z)=mu;
end