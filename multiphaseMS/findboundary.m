function [xdir,ydir,smallu,linind,slinind,m,n,xsizes,ysizes]=findboundary(u,r,disk,M,N,mexed)
[row, column] = find(u);
if isempty(row)
    xdir=0;
    ydir=0;
    smallu=0;
    linind=0;
    slinind=0;
    m=-1;
    n=-1;
    xsizes=0;
    ysizes=0;
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

if mexed
    [xdir,ydir,xsizes,ysizes]=makerowcolmapsz_mex(mask,rl,cl,m,n);
else
    [xdir,ydir]=makerowcolmaps(mask,rl,cl,m,n);
    xsizes=0;
    ysizes=0;
end

[row,col]=find(mask);
linind=sub2ind([M,N],row+miny-1,col+minx-1);
slinind=sub2ind([m,n],row,col);
