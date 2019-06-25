function dist=a2asameparentmetric(alpha1,alpha2)
    Pm=getsymmetries('hex');
    Um=getalpha2alpha();
    numS=size(Pm,3);
    numU=size(Um,3);
    N=size(alpha1,2);
    dist=1.5*ones(1,N);
    
    
    for j=1:numS
        for k=1:numS
            for i=1:numU
                val=abs(alpha2'*Pm(:,:,j)*Um(:,:,i)'*Pm(:,:,k)*alpha1);
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
end

