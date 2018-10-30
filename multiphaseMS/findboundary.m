function [xdir,ydir,smallu,linind,slinind,m,n]=findboundary(u,r,disk,M,N)
[row, column] = find(u);
if isempty(row)
    xdir=0;
    ydir=0;
    smallu=0;
    linind=0;
    slinind=0;
    m=-1;
    n=-1;
    return
end
maxr=max(row);
minr=min(row);
minc=min(column);
maxc=max(column);
miny=max(minr-r,1);
maxy=min(maxr+r,M);
minx=max(minc-r,1);
maxx=min(maxc+r,N);

indy=miny:maxy;
indx=minx:maxx;
smallu=u(indy,indx)*1;
%miny
[m,n]=size(smallu);
mask=findboundarypixels(smallu);
% if sum(mask(:)~=smallu(:))<10
%     xdir=0;
%     ydir=0;
%     smallu=0;
%     linind=0;
%     slinind=0;
%     m=-1;
%     n=-1;
%     return
% end
mask = imdilate(mask,disk);
rl=sum(sum(mask,2)>0);
cl=sum(sum(mask,1)>0);

xdir = cell(1,rl);
step=1;
for i=1:m
    [arrays,count]=decompose(mask(i,:));
    if count>0
        xdir{step}=[i,count,arrays];
        step=step+1;
    end
end
ydir = cell(1,cl);
step=1;
for j=1:n
    [arrays,count]=decompose(mask(:,j));
    if count>0
        ydir{step}=[j,count,arrays];
        step=step+1;
    end
end
[row,col]=find(mask);
linind=sub2ind([M,N],row+miny-1,col+minx-1);
slinind=sub2ind([m,n],row,col);
