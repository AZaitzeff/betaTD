function Rmat=bestparentmatrix(m,n,flatquat,beta,propbeta,onemean)
Rmat=zeros(3,3,m,n);
for i=1:m
    for j=1:n
        ind = sub2ind([m n], i, j);
        if onemean
            [best,~]=bestparentRmat(flatquat(:,ind),propbeta,beta(ind));
        else
            [best,~]=bestparentRmat(flatquat(:,ind),propbeta(:,:,i,j),beta(ind));
        end
        Rmat(:,:,i,j)=best;
    end
end