names={'sim','AFone','RX','AFbig','mapcenter','mapedge','hardbot','hardmid','hardtop','AFbeta'};
addpath('../anglelib/')
full=1;
for i=5:5
name=names{i};
EBSDtemp=load(['../data/' name 'EBSD']);
[M,N]=size(EBSDtemp.CI);
if ~full
    midm=ceil(M/2);
    midn=ceil(N/2);
    lm=max(midm-250,1);
    um=min(midm+250,M);
    ln=max(midn-250,1);
    un=min(midn+250,N);
    suby=lm:um;
    subx=ln:un;
    M=um-lm+1;
    N=un-ln+1;
    EBSD=EBSDtemp.EBSD(suby,subx,:);
    CI=EBSDtemp.CI(suby,subx);
    beta=logical(EBSDtemp.betas(suby,subx));
else
    EBSD=EBSDtemp.EBSD;
    CI=EBSDtemp.CI;
    beta=logical(EBSDtemp.betas);
end
vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
codegen initializeEBSDfast -args {vectorType2,vectorType1,vectorType1l,1,1}
nr=ceil(M/50);
nc=ceil(N/50);
[mapall,dict,kappa,~]=initializeEBSDfast_mex(EBSD,CI,beta,nr,nc);
truebetaEBSD=converttobetamap(EBSD,beta,dict,mapall);
save(['results/' name 'patchfull'],'mapall','dict','kappa','truebetaEBSD');
end