function dist=b2bmetric(beta1,beta2)
    Pm=real(getsymmetries('cubic'));
    numS=size(Pm,3);
    n=size(beta1,1);
    dist=1.5*ones(1,n);
    beta2=real(beta2);
    beta1=real(beta1);
    for j=1:numS
        dist(:,:)=min(dist(:,:),acos(beta2*(Pm(:,:,j)'*beta1')));
    end
    
end