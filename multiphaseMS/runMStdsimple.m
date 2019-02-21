function runMStdsimple(filename,filesave,numpar,num,dt,fid,gs,small)
%fac=40
%gs=50
addpath('../anglelib/')


if nargin<3
    numpar=10;
    num=10;
end


EBSDtemp=load(['../data/' filename 'EBSD.mat']);



addpath('../anglelib/')
dx=1/100;
dy=dx*EBSDtemp.scale;

if small
    [M,N]=size(EBSDtemp.CI);
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
    [M,N]=size(CI);
    beta=logical(EBSDtemp.betas);
end
codegenzaitzeff(M,N);
nr=ceil(M/gs);
nc=ceil(N/gs);
name=['results/' filesave num2str(round(fid))];

timings=zeros(1,num);
if numpar>1
    parpool([1 numpar])
    
    parfor pari=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
        timings(pari)=toc;
    end
    [I,conval,conmap]=confidencemapminmem(name,M,N,num,numpar);
    
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
        timings(i)=toc;
    end
    [I,conval,conmap]=confidencemapminmem(name,M,N,num,numpar);
end
w=5;
[~,bndconval,bndconmap]=probmetric(name,w,num);
vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;


betaEBSD=converttobetamap(EBSD,beta,dict,mapall);
if small
    finalname=['results/' filesave num2str(round(fid)) 'sm'];
else
    finalname=['results/' filesave num2str(round(fid))];
end
save(finalname,'mapall','betaEBSD','dict','energy','conval','conmap','timings','bndconval','bndconmap');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end