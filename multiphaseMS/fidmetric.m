function [val,fidmap]=fidmetric(EBSD,CI,beta,mapall,dict)
[M,N,z]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,z]);
fidmap=zeros(M,N);
EBSDflat=E313toq(EBSDflat);
[K,~]=size(dict);
[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
for k =1:K
    val=0;
    linind=coords(k,1:sizecoords(k));
    bmask=beta(linind);
    alphacoord=linind(~bmask);
    betacoord=linind(bmask);
    if ~isempty(alphacoord)
        fidmap(alphacoord)=360/pi*CI(alphacoord).*(alpbmetric(EBSDflat(alphacoord,:),dict(k,:)));
        val=val+sum(fidmap(alphacoord));
    end
    if ~isempty(betacoord)
        fidmap(betacoord)=360/pi*CI(betacoord).*(b2bmetric(EBSDflat(betacoord,:),dict(k,:)));
        val=val+sum(fidmap(betacoord));
    end
end
val=(val/sum(CI(:)));
