function runMStdpilchank(filename,filesave,dt,fids,small,numwanted)
%fac=40
%gs=50



num=numel(fids);
EBSDtemp=load(['data/' filename 'EBSD.mat']);


dx=1/100;
dy=dx*EBSDtemp.scale;

if small
    [M,N]=size(EBSDtemp.CI);
    midm=ceil(M/2);
    midn=ceil(N/2);
    lm=max(midm-small,1);
    um=min(midm+small,M);
    ln=max(midn-small,1);
    un=min(midn+small,N);
    suby=lm:um;
    subx=ln:un;
    M=um-lm+1;
    N=un-ln+1;
    EBSD=EBSDtemp.EBSD(suby,subx,:);
    CI=EBSDtemp.CI(suby,subx);
    IQ=EBSDtemp.IQ(suby,subx);
    beta=logical(EBSDtemp.betas(suby,subx));
else
    EBSD=EBSDtemp.EBSD;
    CI=EBSDtemp.CI;
    IQ=EBSDtemp.IQ;
    [M,N]=size(CI);
    beta=logical(EBSDtemp.betas);
end
 CI=CIfunc(CI);

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
codegen zpilchank -args {vectorType2,vectorType1,vectorType1,vectorType1l,1,1,1,1,1,1}
tic;
[newmap,newbetadict,K]=zpilchank_mex(EBSD,IQ,CI,beta,3,20,3,3,1,numwanted);
timing=toc;
save(['results/' filesave 'P'],'newmap','newbetadict','timing');
kappa=ones(1,K)*50;
dtstop=2^(-10);
between=2;
nt=6;

vectorType3 = coder.typeof(1, [4 K], [false true]);
vectorType4 = coder.typeof(1, [1 K], [false true]);
codegen EBSDimgseggivenint -args {vectorType1,vectorType3,vectorType4,vectorType2,vectorType1,vectorType1l,1,1,1,1,1}
for i=1:num
    fid=fids(i);
    tic;
    [mapall,dict,energy,gsizes,flag]=EBSDimgseggivenint_mex(newmap,newbetadict,kappa,EBSD,CI,beta,fid,dt,dx,dy,400);
    timing=toc;
    finalname=['results/' filesave num2str(round(fid))];
    betaEBSD=converttobetamap(EBSD,beta,dict,mapall);
    save(finalname,'mapall','betaEBSD','dict','energy','gsizes','fid','timing','flag');
end
