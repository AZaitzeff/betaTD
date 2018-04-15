function run2pMS(filename,filesave,fids,biggest,numpar,betathres)
EBSDtemp=load(['../data/' filename 'EBSD.mat']);
[x,y]=size(EBSDtemp.CI);
d=ceil(x/biggest);
xranges=ones(1,d+1);
for i=1:d
    xranges(i+1)=ceil(x/d*(i));
end
d=ceil(y/biggest);
yranges=ones(1,d+1);
for i=1:d
    yranges(i+1)=ceil(y/d*(i));
end
m=numel(xranges)-1;
n=numel(yranges)-1;
total=m*n;
parpool(numpar)
factor=100;
parfor section=1:total
    for fid=fids
        EBSDtemp=load(['../data/' filename 'EBSD.mat']);
        addpath('../anglelib/')
        mapall=load(['../data/' filename 'ff.mat']);
        i=mod(section-1,m)+1;
        j=ceil(section/m);
        rows=(max(xranges(i)-factor,1):min(xranges(i+1)+factor,xranges(end)));
        cols=(max(yranges(j)-factor,1):min(yranges(j+1)+factor,yranges(end)));

        EBSD=EBSDtemp.EBSD(rows,cols,:);
        CI=EBSDtemp.CI(rows,cols);
        betas=EBSDtemp.betas(rows,cols);
        map=mapall.map(rows,cols);
        
        tic;
        twophaseMS(EBSD,CI,betas,map,fid,betathres,section,filename);
        toc;

    end
end
poolobj = gcp('nocreate');
delete(poolobj);
puttogether(filename,filesave,fids,betathres,xranges,yranges,factor)
end