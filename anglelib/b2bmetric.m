function dist=b2bmetric(beta1,beta2)
    Pm=getsymmetries('cubic');
    numS=size(Pm,3);
    N=size(beta1,1);
    n=size(beta2,1);
    dist=1.5*ones(N,n);
    
    for j=1:numS
        dist(:,:)=min(dist(:,:),acos(abs((beta1')'*(Pm(:,:,j)*beta2')))');
    end
end