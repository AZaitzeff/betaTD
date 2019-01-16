function dist=b2bmetric(beta1,beta2)
    Pm=getsymmetries('cubic');
    numS=size(Pm,3);
    n=size(beta1,1);
    dist=1.5*ones(1,n);

    for j=1:numS
        for i=1:N
            dist(i)=min(dist(i),acos(abs((beta2*Pm(:,:,j)')*beta1(i,:)')));
        end
    end
    
end