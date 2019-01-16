function dist=b2bmetric(beta1,beta2)
    Pm=getsymmetries('cubic');
    numS=size(Pm,3);
    n=size(beta1,1);
    dist=1.5*ones(1,n);

    for j=1:numS
        val=abs((beta2*Pm(:,:,j)')*beta1');
        if any(val>1)
            val(val>1)-1
        end
        dist=min(dist,acos(abs((beta2*Pm(:,:,j)')*beta1')));
    end
    
end