function dist=b2bmetric(beta1,beta2)
    Pm=getsymmetries('cubic');
    numS=size(Pm,3);
    n=size(beta1,1);
    dist=1.5*ones(1,n);

    for j=1:numS
        if ~isreal((beta2*Pm(:,:,j)')*beta1')
            (beta2*Pm(:,:,j)')*beta1'
        end
        if ~isreal(abs((beta2*Pm(:,:,j)')*beta1'))
            abs((beta2*Pm(:,:,j)')*beta1')
        end
        dist=min(dist,acos(abs((beta2*Pm(:,:,j)')*beta1')));
    end
    
end