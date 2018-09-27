function energy=EBSDtdE(u,EBSD,CI,dict,fid)
energy=0;
dt=.01;


K=size(u,2);
[m,n]=size(u{1});%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[m*n,1]);


[ls]=td2dz(u,dt);

for k=1:K
    ind=(u{k}(:)>0);
    if sum(ind)>0
        energy=energy+2/sqrt(dt)*sum(ls{k}(~ind));
        energy=energy+sum(fid*CIflat(ind(:)).*alpbmetric(EBSDflat(ind(:),:),dict{k})');
    end
end



end