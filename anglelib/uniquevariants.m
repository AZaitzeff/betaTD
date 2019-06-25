function [num]=uniquevariants(beta,alpha,tol)
    T=betatoalphatrans();
    numB=size(T,3);
    num=zeros(1,numB);
    for i=1:numB
        candalpha=(beta'*T(:,:,i))';
        num(i)=sum(a2ametric(alpha,candalpha)<tol);
    end
end