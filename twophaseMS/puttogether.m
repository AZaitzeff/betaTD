function puttogether(filename,fids,xranges,yranges)
addpath('../anglelib/')
numrows=xranges(end);
numcols=yranges(end);
m=numel(xranges);
n=numel(yranges);
total=m*n;
for fid=fids
    
mapall=zeros(numrows,numcols);
betaEBSD=zeros(numrows*numcols,4);
load(['../data/' filename 'alpha.mat']);
z=1;
for section=1:total
    
    load(['results/' filename 'part' num2str(section) num2str(fid)]);
    
    factor=25;


    i=mod(section-1,m)+1;
    j=ceil(section/m);
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
                 
    mapall(rowssmall,colssmall)=temp(factor+1-(i==1)*factor:factor+150-(i==1)*factor,factor+1-(j==1)*factor:factor+175-(j==1)*factor);    

end
[m,n,z]=size(EBSD);
EBSDflat=reshape(EBSD, [m*n,z]);
EBSDflat=E313toq(EBSDflat);
CIflat=reshape(CI, [m*n,1]);
[clusterlist,~,~] = unique(mapall);
T=alphatobetatrans();
Pm=getsymmetries('cubic');
Pall=zeros(4,4,144);
for i=1:144
    Pall(:,:,i)=T(:,:,mod(i-1,6)+1)'*Pm(:,:,ceil(i/6));
end
for z=clusterlist'
    zfun=find(z==mapall);
    [mu, ~, ~, ~] = VMFEM(EBSDflat(zfun,:), Pall,CIflat(zfun),1,16);
    for i=1:4
        betaEBSD(zfun,i)=mu(i);
    end
end
betaEBSD=qtoE313(betaEBSD);
betaEBSD=reshape(betaEBSD,[m,n,3]);
save(['results/' filename num2str(fid)],'mapall','betaEBSD');
end
end
    
