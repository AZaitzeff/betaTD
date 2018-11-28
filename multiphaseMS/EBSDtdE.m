function energy=EBSDtdE(mapall,EBSD,CI,dict,fid,dt,dx,dy,enec)
energy=0;

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


if enec
     ene=zeros(K,K);
     for k1=1:K
         for k2=k1+1:K
             energy=readshockley(dict(k1,:),dict(k2,:));
             
             ene(k1,k2)=energy;
             ene(k2,k1)=energy;
         end
     end
     

     
    for k=1:K
        ind=(u{k}(:)>0);
        energy=energy+sum(fid*CIflat(ind(:)).*alpbmetric(EBSDflat(ind(:),:),dict(k,:))');
        if sum(ind)>0
            for k2=1:K
                if k2~=k
                    energy=energy+ene(k,k2)*1/sqrt(dt)*sum(ls{k2}(ind));
                end
            end
        end
     end
else
    
for k=1:K
    ind=(u{k}(:)>0);
    if sum(ind)>0
        energy=energy+2/sqrt(dt)*sum(ls{k}(~ind));
        energy=energy+sum(fid*CIflat(ind(:)).*alpbmetric(EBSDflat(ind(:),:),dict(k,:))');
    end
end 
     
end






end