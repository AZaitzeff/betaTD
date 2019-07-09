function [linind,slinind,fullsz,bnds]=findboundary(w,bnds,xbc,ybc,M,N)
minr=bnds(1);
maxr=bnds(2);
minc=bnds(3);
maxc=bnds(4);
miny=max(minr-(w+2),1);
maxy=min(maxr+(w+2),M);
minx=max(minc-(w+2),1);
maxx=min(maxc+(w+2),N);
m=maxy-miny+1;
n=maxx-minx+1;
%indy=miny:maxy;
%indx=minx:maxx;
bnds(1)=miny;
bnds(2)=maxy;
bnds(3)=minx;
bnds(4)=maxx;
%smallu=u(indy,indx)==k;
%SE = strel('disk',w,4);
%J = imdilate(smallu,SE);
%[row,col]=find(J);
%W=smallu*0-1;
x=xbc-minx+1;
y=ybc-miny+1;

[row,col,fullsz]=zgrow(y,x,w,m,n);
%[row,col] = pgrow3(int32(y),int32(x),w,int32(W));
linind=sub2ind([M,N],row+miny-1,col+minx-1);
slinind=sub2ind([m,n],row,col);



