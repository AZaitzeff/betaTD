function smoothed=smoothonegrain(EBSD,CI,beta)


[m,n,~]=size(EBSD);
smallflatq=E313toq(reshape(EBSD,m*n,3)');
[mu, ~] = estimatebeta(m*n, CI(:), smallflatq, beta(:), [], [], 400);
propbeta = RMatOfQuat(mu);
Rmat = bestparentmatrix(m,n,smallflatq,beta,propbeta,1);

uin=zeros(3,3,m,n);
for i=1:m
    for j=1:n
        uin(:,:,i,j)=propbeta;
    end
end
%your code here
[smoothed,~,~]=so3implicitfid(100,1/100,Rmat,uin,100);
end

