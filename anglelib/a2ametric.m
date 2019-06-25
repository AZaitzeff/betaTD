function dist=a2ametric(alpha1,alpha2)
    Pm=getsymmetries('hex');
    numS=size(Pm,3);
    N=size(alpha1,2);
    dist=1.5*ones(1,N);
    
    for j=1:numS
    	val=abs(alpha2'*Pm(:,:,j)*alpha1);
        for n=1:N
            if val(n)>=1
                dist(n)=0;
            else
                dist(n)=min(dist(n),2*acos(val(n)));
            end
        end
    end
end