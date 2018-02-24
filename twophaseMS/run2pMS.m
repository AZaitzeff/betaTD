function run2pMS(filename,fids,factor)
EBSDtemp=load(['../data/' filename 'alpha.mat']);
[x,y]=size(EBSDtemp.CI);
d=ceil(x/factor);
xranges=ones(1,d+1);
for i=1:d
    xranges(i+1)=ceil(x/d*(i));
end
d=ceil(y/factor);
yranges=ones(1,d+1);
for i=1:d
    yranges(i+1)=ceil(y/d*(i));
end
m=numel(xranges)-1;
n=numel(yranges)-1;
total=m*n;
parpool(10)
%fids=[75,100,125,150,175];
parfor section=1:total
    for fid=fids
        EBSDtemp=load(['../data/' filename 'alpha.mat']);
        addpath('../anglelib/')
        mapall=load(['../data/' filename 'ff.mat']);
        factor=25;
        i=mod(section-1,m)+1;
        j=ceil(section/m);
        rows=(max(xranges(i)-factor,1):min(xranges(i+1)+factor,xranges(end)));
        cols=(max(yranges(j)-factor,1):min(yranges(j+1)+factor,yranges(end)));

        EBSD=EBSDtemp.EBSD(rows,cols,:);
        CI=EBSDtemp.CI(rows,cols);
        map=mapall.map(rows,cols);
        tic;
        twophaseMS(EBSD,CI,map,fid,section,filename);
        toc;

    end
end

puttogether(filename,fids,xranges,yranges)
end