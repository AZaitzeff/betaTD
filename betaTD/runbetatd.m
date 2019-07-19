function runbetatd(filename,filesave,numpar,num,fid,dt,gs)
%fac=40
%gs=50

if nargin<3
    numpar=10;
    num=10;
end


EBSDtemp=load(['data/' filename 'EBSD.mat']);




EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);

CI=CIfunc(CI);
codegenzaitzeff(M,N);
nr=ceil(M/gs);
nc=ceil(N/gs);
name=['results/' filesave num2str(round(fid))];

timings=zeros(1,num);
scale=EBSDtemp.scale;
if numpar>1
    parpool([1 numpar])
    
    parfor pari=1:num
        tic;
        dttemp=dt*2^(rand()-1/2);
        MStd(EBSD,CI,beta,fid,filesave,dttemp,scale,nr,nc,pari);
        timings(pari)=toc;
    end
    
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        dttemp=dt*2^(rand()-1/2);
        MStd(EBSD,CI,beta,fid,filesave,dttemp,scale,nr,nc,i);
        timings(i)=toc;
    end
end
w=1;
[I,conval,conmap]=confidencemapminC(name,M,N,num);
[~,bndconval,bndconmap]=probmetric(name,w,num);
vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;


betaEBSD=converttobetamap(EBSD,beta,dict,mapall);

    
finalname=['results/' filesave num2str(round(fid))];

save(finalname,'mapall','betaEBSD','dict','energy','conval','conmap','timings','bndconval','bndconmap');

for i=1:num
 delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end