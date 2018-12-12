function runMStd(filename,filesave,numpar,num,smdiam)
addpath('../anglelib/')

fid=200;
dt=2^-6;
if nargin<3
    numpar=10;
    num=10;
end
if nargin<5
    smdiam=50;
end


EBSDtemp=load(['../data/' filename 'EBSD.mat']);
addpath('../anglelib/')
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
beta=logical(EBSDtemp.betas);
codegenzaitzeff(M,N);


Ks=[ceil(N/smdiam*1.2),ceil(M/smdiam*1.2)];
dx=1/(2*smdiam);
dy=dx*EBSDtemp.scale;

%betas=EBSDtemp.betas(rows,cols);
timings=zeros(1,num);

name=['results/' filesave num2str(round(fid))];
smallK=ceil((Ks(1)*Ks(2))/8);

if numpar>1
    parpool(numpar)
    
    parfor i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,Ks,filesave,dt,dx,dy,i);
        timings(i)=toc;
    end
    pause(1);
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
    
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
else
    for i=1:num
        tic;
        MStd(EBSD,CI,beta,fid,Ks,filesave,dt,dx,dy,i);
        timings(i)=toc;
    end
    pause(1);
    [I,conval,conmap]=confidencemap(name,M,N,smallK,num,numpar);
end


vars=load(['results/' filesave num2str(round(fid)) num2str(I)]);



mapall=vars.mapall;
dict=vars.dict;
energy=vars.energy;


betaEBSD=zeros(M*N,4);
K = size(dict,1);
alpha=reshape(EBSD, [M*N,3]);
alpha=E313toq(alpha);

[~,~,~,coords,sizecoords,~]=  bndcoords(mapall,K);
for k=1:K
    mu=dict(k,:);
    indices=coords(k,1:sizecoords(k));
    mask=beta(indices);
    acoord=indices(~mask);
    [best,~]=alpbbeta(alpha(acoord,:),mu);
    betaEBSD(acoord,:)=best;
    bcoord=indices(mask);
    betaEBSD(bcoord,:)=alpha(bcoord,:);

    
end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[M,N,3]);
save(['results/' filesave num2str(round(fid))],'mapall','betaEBSD','dict','energy','conval','conmap','timings');

for i=1:num
    delete(['results/' filesave num2str(round(fid)) num2str(i) '.mat']);
end

