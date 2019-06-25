function dist=alpbmetric(alpha,beta)
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    numB=size(T,3);
    numS=size(Pm,3);
    N=size(alpha,2);
    dist=1.5*ones(1,N);
    
    for i=1:numB
        for j=1:numS
            val=abs((beta'*Pm(:,:,j))*(alpha'*T(:,:,i))');
            for n=1:N
                if val(n)>=1
                    dist(n)=0;
                else
                    dist(n)=min(dist(n),2*acos(val(n)));
                end
            end
        end
    end
end