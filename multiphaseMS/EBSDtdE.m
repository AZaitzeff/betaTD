function energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dx,dy)
energy=0;
dt=1/2^12;

K=max(mapall(:));

[m,n]=size(mapall);%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[m*n,1]);
u=cell(1,K);
    
for k=1:K
    u{k}=mapall==k;
end

[ls]=td2dz(u,dt,dx,dy);

for k=1:K
    ind=(u{k}(:)>0);
    if sum(ind)>0
        energy=energy+2/sqrt(dt)*sum(ls{k}(~ind));
        energy=energy+sum(fid*CIflat(ind(:)).*alpbmetric(EBSDflat(ind(:),:),dict{k})');
    end
end



end