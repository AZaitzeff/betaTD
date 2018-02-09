function dist=abmetric2(alpha,beta)
    T=betatoalphatrans();
    Pm=getsymmetries('hex');
    numB=size(T,3);
    numS=size(Pm,3);
    N=size(beta,1);
    n=size(alpha,1);
    dist=1.5*ones(N,n);
    
    for i=1:numB
        for j=1:numS
            dist(:,:)=min(dist(:,:),acos(abs((T(:,:,i)*beta')'*(Pm(:,:,j)*alpha'))));
        end
    end
end