function [newg1,kap]=estimatebeta(N,CI,EBSDflat,beta,g1,kappa,numsub)
    T=alphatobetatrans();%load the symmetries
    Pm=getsymmetries('cubic');
    Pall=zeros(4,4,144);
    for i=1:144
        Pall(:,:,i)=Pm(:,:,ceil(i/6))*T(:,:,mod(i-1,6)+1)';
    end
    newind=datasamplez(1:N,numsub,CI);
    EBSDtemp=EBSDflat(:,newind);
    cmask=beta(newind);
    alphaEBSD=EBSDtemp(:,~cmask);
    betaEBSD=EBSDtemp(:,cmask);
    [newg1, kap, ~] = VMFEMzfast(alphaEBSD, Pall,betaEBSD, Pm,1,g1,kappa);
    vala=alpbmetric(alphaEBSD,newg1);
    valb=b2bmetric(betaEBSD,newg1);
    tol=median([vala valb]);
    maska=vala<tol;
    maskb=valb<tol;
    alpha2=alphaEBSD(:,maska);
    beta2=betaEBSD(:,maskb);
    [newg1, kap, ~] = VMFEMzfast(alpha2, Pall,beta2, Pm,1,newg1,kap);
end