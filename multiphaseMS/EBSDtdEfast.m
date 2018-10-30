function energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dx,dy)
energy=0;
dt=1/2^12;
K=max(mapall(:));
[M,N]=size(mapall);%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI,[M*N,1]);

for k=1:K
    ind=(mapall(:)==k);
    if sum(ind)>0
        energy=energy+sum(fid*CIflat(ind).*alpbmetric(EBSDflat(ind,:),dict{k})');
    end
end

r=ceil(200*sqrt(dt));
disk=strel('disk',r-1,4);
for k=1:K
    temp=(mapall==k);
    [xdir,ydir,smallu,~,~,hor,ver]=findboundary(temp,r,disk,M,N);
    newu=ADI(smallu,dt,dx,dy,xdir,ydir,hor,ver);
    mask=(smallu(:)==0);
    energy=energy+2/sqrt(dt)*sum(newu(mask));
end



end