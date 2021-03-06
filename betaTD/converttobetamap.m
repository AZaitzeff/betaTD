function betaEBSD=converttobetamap(EBSD,beta,dict,mapall)
[M,N]=size(mapall);
betaEBSD=zeros(M*N,4);
K = size(dict,2);
alpha=reshape(EBSD, [M*N,3]);
alpha=E313toq(alpha');

[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
for k=1:K
    mu=dict(:,k);
    indices=coords(k,1:sizecoords(k));
    mask=beta(indices);
    acoord=indices(~mask);
    [best,~]=alpbbeta(alpha(:,acoord),mu);
    betaEBSD(acoord,:)=best;
    bcoord=indices(mask);
    betaEBSD(bcoord,:)=alpha(:,bcoord)';

    
end

betaEBSD=qtoE313(betaEBSD')';
betaEBSD=reshape(betaEBSD,[M,N,3]);