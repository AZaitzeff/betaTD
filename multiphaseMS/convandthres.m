function [newmapall,curmin]=convandthres(mapall,curmin,dict,CI,EBSD,K,active,dx,dy,dt,fid,rmspots,mexed)
if nargin<12
    rmspots=0;
end
if nargin<13
    mexed=0;
end

[M,N]=size(mapall);
r=ceil(200*sqrt(dt));
disk=strel('disk',r-1,4);
newmapall=mapall;
for k=1:K
    if active(k)
        temp=(mapall==k);
        if rmspots>1
            temp=bwareaopen(temp,rmspots);
        end
        [xdir,ydir,smallu,linind,slinind,m,n,xsizes,ysizes]=findboundary(temp,r,disk,M,N,mexed);
        if m>0
            if mexed
                newu=ADIz_mex(smallu,dt,dx,dy,xdir,ydir,xsizes,ysizes,m,n);
            else
                newu=ADI(smallu,dt,dx,dy,xdir,ydir,m,n);
            end
            S=2/sqrt(dt)*(-newu(slinind))+fid*CI(linind).*alpbmetric(EBSD(linind,:),dict{k})';
            [minval,I]=findminz(curmin(linind),S);
            curmin(linind(I))=minval(I);
            newmapall(linind(I))=k;
        end
    end
end