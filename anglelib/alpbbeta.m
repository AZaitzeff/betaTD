function [best,dist]=alpbbeta(alpha,beta)
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    numB=size(T,3);
    numS=size(Pm,3);
    n=size(alpha,1);
    dist=1.5*ones(1,n);
    best=zeros(n,4);
    
    for i=1:numB
        for j=1:numS
            tempdist=acos(abs((Pm(:,:,j)'*T(:,:,i)*alpha')'*(beta')))';
            tempbetas=(Pm(:,:,j)'*T(:,:,i)*alpha')';
            best(tempdist<dist,:)=tempbetas(tempdist<dist,:);
            dist(tempdist<dist)=tempdist(tempdist<dist);
        end
    end
end