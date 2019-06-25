function [betasol,unique]=bestbeta(alpha,beta,tol)
    if tol>.75
        tol=tol/180*pi;
    end
    T=alphatobetatrans();
    numB=size(T,3);
    numberbeta=zeros(1,numB);
    thebetas=zeros(4,numB);
    mainalpha=alpha(:,1);
    for i=1:numB
        candbeta=(mainalpha'*T(:,:,i))';
        thebetas(:,i)=candbeta;
        numberbeta(i)=(sum(alpbmetric(alpha,candbeta)<tol)+sum(b2bmetric(beta,candbeta)<tol));
    end
    numberbeta=numberbeta/max(numberbeta);
    
    [WhichBetaGrain] = find(numberbeta == 1);
if numel(WhichBetaGrain) == 1  % there is only one solution!
    betasol = thebetas(:,WhichBetaGrain);
    unique = 1;
else
    betasol = thebetas(:,WhichBetaGrain);
    unique = 0;
end
end