function runMStdsimple(filename,filesave,numpar,num,dt,fid,gs,fac)
%fac=40
%gs=50
addpath('../anglelib/')


if nargin<3
    numpar=10;
    num=10;
end


EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
dx=1/(2*gs);
dy=dx*EBSDtemp.scale;

EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N,0);
nr=ceil(M/fac*gs/50);
nc=ceil(N/fac*gs/50);
name=['results/' filesave num2str(round(fid))];

timings=zeros(1,num);
smallK=ceil((nr*nc)/8);
if numpar>1
    parpool(numpar)
    
    parfor pari=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,pari);
        timings(pari)=toc;
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,filesave,dt,dx,dy,nr,nc,i);
        timings(i)=toc;
    end
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
end

vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;


betaEBSD=converttobetamap(EBSD,beta,dict,mapall);
save(['results/' filesave num2str(round(fid))],'mapall','betaEBSD','dict','energy','conval','conmap','timings');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end