function [xdir,ydir,xsizes,ysizes,smallu,linind,slinind,bndcoor,bndsz,fullsz,m,n]=findboundary(u,k,w,bnds,xbc,ybc,M,N)
minr=bnds(1);
maxr=bnds(2);
minc=bnds(3);
maxc=bnds(4);

miny=max(minr-(w+5),1);
maxy=min(maxr+(w+5),M);
minx=max(minc-(w+5),1);
maxx=min(maxc+(w+5),N);
m=maxy-miny+1;
n=maxx-minx+1;
indy=miny:maxy;
indx=minx:maxx;
smallu=u(indy,indx)==k;
%SE = strel('disk',w,4);
%J = imdilate(smallu,SE);
%[row,col]=find(J);
%W=smallu*0-1;
x=xbc-minx+1;
y=ybc-miny+1;

[row,col,ind]=zgrow(y,x,w,m,n);
%[row,col] = pgrow3(int32(y),int32(x),w,int32(W));
linind=sub2ind([M,N],row+miny-1,col+minx-1);
bndsz=ind(1);
fullsz=ind(2);
bndcoor=sort(linind(1:bndsz));
slinind=sub2ind([m,n],row,col);
mask=smallu*0;
mask(slinind)=1;
[xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,m,n);


