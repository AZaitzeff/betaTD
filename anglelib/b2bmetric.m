function dist=b2bmetric(beta1,beta2)
    Pm=getsymmetries('cubic');
    numS=size(Pm,3);
    N=size(beta1,1);
    dist=1.5*ones(1,N);

    for j=1:numS
        val=abs((beta2*Pm(:,:,j)')*beta1');
        for n=1:N
            if val(n)>=1
                dist(n)=0;
            else
                dist(n)=min(dist(n),acos(val(n)));
            end
        end
    end
    
end