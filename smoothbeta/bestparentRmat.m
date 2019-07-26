function [best,dist]=bestparentRmat(alpha,beta,isbeta)
    alpha=QuatOfRMat(alpha);
    beta=QuatOfRMat(beta);
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    numB=size(T,3);
    numS=size(Pm,3);
    dist=1.5;
    best=zeros(1,4);
    if isbeta
        for j=1:numS
            tempdist=2*acos(abs((beta'*Pm(:,:,j))*(alpha)));
            tempbetas=alpha'*Pm(:,:,j)';
            if tempdist<dist
                best=tempbetas;
                dist=tempdist;
            end
        end
        
    else
        for i=1:numB
            for j=1:numS
                tempdist=2*acos(abs((beta'*Pm(:,:,j))*(alpha'*T(:,:,i))'));
                tempbetas=alpha'*T(:,:,i)*Pm(:,:,j)';
                if tempdist<dist
                    best=tempbetas;
                    dist=tempdist;
                end
            end
        end
    end
    best=RMatOfQuat(best);
    
end