smallflatq=E313toq(reshape(EBSD,300*150,3)');

[m,n,~]=size(EBSD);
[mu,kap]=estimatebeta(m*n,CI(:),EBSD,beta(:),[],[],400);
meanbeta=RMatOfQuat(mu);
Rmat=zeros(3,3,m,n);
for i=1:m
    for j=1:n
        ind = sub2ind([m n], i, j);
        [best,dist]=bestparentRmat(smallflatq(:,ind),meanbeta,beta(ind));
        Rmat(:,:,i,j)=best;
    end
end


%your code here