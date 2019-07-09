function [vals2]=matchmetric2(EBSD, CI,beta,mapall,dict)

[K,~]=size(dict);
[neighbors]=  findneigh(mapall,K);
vals2=zeros(1,0);
z=1;

[M,N,ele]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,ele]);
EBSDflat=E313toq(EBSDflat);
[K,~]=size(dict);

[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
for k =1:K
    val=0;
    linind=coords(k,1:sizecoords(k));
    bmask=beta(linind);
    alphacoord=linind(~bmask);
    betacoord=linind(bmask);
    for j=1:K
        if neighbors(k,j)
            if ~isempty(alphacoord)
                val=val+sum(CI(alphacoord).*(alpbmetric(EBSDflat(alphacoord,:),dict(j,:))));
            end
            if ~isempty(betacoord)
                val=val+sum(CI(betacoord).*(b2bmetric(EBSDflat(betacoord,:),dict(j,:))));
            end
            vals2(z)=val/sum(CI(linind));
            z=z+1;
        end
    end
end

vals2=vals2*2*180/pi;