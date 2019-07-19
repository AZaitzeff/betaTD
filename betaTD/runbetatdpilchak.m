function runbetatdpilchak(filename,filesave,sizepar,mex)
% sizepar should be:
% 1: grains around 50 by 50 or smaller
% 2: grains around 100 by 100  (a good pick if you do not know)
% 3: grains around 200 by 200 or larger

EBSDtemp=load(['data/' filename 'EBSD.mat']);
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
IQ=EBSDtemp.IQ;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
 CI=CIfunc(CI);
 
    if sizepar==1
        fids=[100,150,200,250];
        dt=2^-5;
        numwanted=ceil(20*(M*N)/50^2);
    elseif sizepar==2
        fids=[75,100,125,150,175];
        dt=2^-4;
        numwanted=ceil(20*(M*N)/100^2);
    else
        fids=[25,50,75,100];
        dt=2^-3;
        numwanted=ceil(20*(M*N)/200^2);
    end


num=numel(fids);
scale=EBSDtemp.scale;



if mex
vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
codegen zpilchak -args {vectorType2,vectorType1,vectorType1,vectorType1l,1,1,1,1,1,1}
    tic;
    [newmap,newbetadict,K]=zpilchak_mex(EBSD,IQ,CI,beta,5,10,3,3,1,numwanted);
    timing=toc;
    save(['results/' filesave 'P'],'newmap','newbetadict','timing');
    kappa=ones(1,K)*50;
else
    tic;
    [newmap,newbetadict,K]=zpilchak(EBSD,IQ,CI,beta,5,10,3,3,1,numwanted);
    timing=toc;
    save(['results/' filesave 'P'],'newmap','newbetadict','timing');
    kappa=ones(1,K)*50;

end

if mex
    vectorType3 = coder.typeof(1, [4 K], [false true]);
    vectorType4 = coder.typeof(1, [1 K], [false true]);
    codegen EBSDimgseggivenint -args {vectorType1,vectorType3,vectorType4,vectorType2,vectorType1,vectorType1l,1,1,1}
end
for i=1:num
    fid=fids(i);
    if mex
        tic;
        [mapall,dict,energy,gsizes]=EBSDimgseggivenint_mex(newmap,newbetadict,kappa,EBSD,CI,beta,fid,dt,scale);
        timing=toc;
    else
        tic;
        [mapall,dict,energy,gsizes]=EBSDimgseggivenint(newmap,newbetadict,kappa,EBSD,CI,beta,fid,dt,scale);
        timing=toc;
    end
    finalname=['results/' filesave num2str(round(fid))];
    betaEBSD=converttobetamap(EBSD,beta,dict,mapall);
    save(finalname,'mapall','betaEBSD','dict','energy','gsizes','fid','timing');
end
