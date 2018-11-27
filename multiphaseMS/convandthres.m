function [newmapall,curmin]=convandthres(mapall,curmin,dict,CI,EBSD,K,active,dx,dy,dt,fid,rmspots,mexed)
if nargin<12
    rmspots=0;
end
if nargin<13
    mexed=0;
end

[M,N]=size(mapall);
fac=1/(50*(dx+dy));
w=ceil(fac*600*sqrt(dt));
newmapall=mapall;
[xbdcoor,ybdcoor,sizeboor,minmaxrowcol]=  bndcoords(mapall,K);
for k=1:K
    if active(k)
        total=sizeboor(k);
        [xdir,ydir,smallu,linind,slinind,m,n,xsizes,ysizes]=findboundary(mapall,k,w,minmaxrowcol(k,:),xbdcoor(k,1:total)',ybdcoor(k,1:total)',M,N,rmspots,mexed);
        if m>0
            if mexed
                newu=ADIz_mex(smallu,dt,dx,dy,xdir,ydir,xsizes,ysizes,m,n);
            else
                newu=ADI(smallu,dt,dx,dy,xdir,ydir,m,n);
            end
            mask=(newu(slinind)>.025 & newu(slinind)<.995);
            slinind=slinind(mask);
            linind=linind(mask);
            S=2/sqrt(dt)*(-newu(slinind))+fid*CI(linind).*alpbmetric(EBSD(linind,:),dict{k})';
            [minval,I]=findminz(curmin(linind),S);
            curmin(linind(I))=minval(I);
            newmapall(linind(I))=k;
        end
    end
end