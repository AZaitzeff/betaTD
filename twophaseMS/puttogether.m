function puttogether(filename,fids,xranges,yranges,factor)
addpath('../anglelib/')
numrows=xranges(end);
numcols=yranges(end);
mr=numel(xranges)-1;
nr=numel(yranges)-1;
total=mr*nr;
for fid=fids
    
mapall=zeros(numrows,numcols);
betaEBSD=zeros(numrows*numcols,4);

z=1;
for section=1:total
    
    load(['results/' filename 'part' num2str(section) num2str(fid)]);
    

    i=mod(section-1,mr)+1;
    j=ceil(section/mr);
    rows=(max(xranges(i)-factor,1):min(xranges(i+1)+factor,xranges(end)));
    cols=(max(yranges(j)-factor,1):min(yranges(j+1)+factor,yranges(end)));
    rowssmall=max(xranges(i),1):min(xranges(i+1),xranges(end));
    colssmall=max(yranges(j),1):min(yranges(j+1),yranges(end));
    smallmap=mapall(rows,cols);
    temp=smallmap*0;
    [clusterlist,~,~] = unique(map);
    [clusterlist2,~,~] = unique(smallmap);
    for part1=clusterlist'
        maxval=30;
        argmaxval=0;
        for part2=clusterlist2'
            if part2~=0
                candmaxval=sum(sum((smallmap==part2).*(map==part1)));
                if candmaxval>maxval
                    maxval=candmaxval;
                    argmaxval=part2;
                end
            end
        end
        if argmaxval
            temp(map==part1)=argmaxval;
        else
            temp(map==part1)=z;
            z=z+1;
        end
    end
    mapall(rowssmall,colssmall)=temp(factor+1-(i==1)*factor:factor+xranges(i+1)-xranges(i)+1-(i==1)*factor,factor+1-(j==1)*factor:factor+yranges(j+1)-yranges(j)+1-(j==1)*factor);    

end
EBSDtemp=load(['../data/' filename 'EBSD.mat']);
EBSD=EBSDtemp.EBSD;
[m,n,~]=size(EBSD);
CI=EBSDtemp.CI;
betas=EBSDtemp.betas;
betas(:)=0;
[dict]=estimatebetas(EBSD,CI,betas,mapall,[],0,1,500);
[clusterlist,~,~] = unique(mapall);

for z=clusterlist'
    mu=dict(z);
    zfun=(z==mapall);
    for i=1:4
        betaEBSD(zfun,i)=mu(i);
    end
end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[m,n,3]);
save(['results/' filename num2str(fid)],'mapall','betaEBSD');
end
end
    
