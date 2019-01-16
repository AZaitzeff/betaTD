function dist=alpbmetric(alpha,beta)
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    numB=size(T,3);
    numS=size(Pm,3);
    n=size(alpha,1);
    dist=1.5*ones(1,n);
    
    for i=1:numB
        for j=1:numS
            dist(:,:)=min(dist(:,:),acos(abs((beta*Pm(:,:,j)'*T(:,:,i))*alpha')));
        end
    end
end