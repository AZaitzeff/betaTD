function energy=EBSDtdEfast(mapall,EBSD,CI,dict,fid,dt,dx,dy)
energy=0;
K=max(mapall(:));
[M,N]=size(mapall);%size of level set function
[~,~,z]=size(EBSD);
EBSDflat=reshape(EBSD,[M*N,z]);
EBSDflat=E313toq(EBSDflat);


[xbdcor,ybdcor,sizebdcor,coords,sizecor,minmaxrowcol]=  bndcoords(mapall,K);
for k=1:K
    linind=coords(k,1:sizecor(k));
    energy=energy+sum(fid*CI(linind).*alpbmetric(EBSDflat(linind,:),dict(k,:)));  
end
w=ceil(fac*600*sqrt(dt));
for k=1:K
    

        total=sizebdcor(k);
        [xdir,ydir,xsizes,ysizes,smallu,~,~,m,n]=findboundary(newmapall,k,w,minmaxrowcol(k,:),xbdcor(k,1:total)',ybdcor(k,1:total)',M,N);



        newu=TSz(smallu*1,dt,10,dx,dy,xdir,ydir,xsizes,ysizes,m,n);
        %mask=(newu(slinind)>.025 & newu(slinind)<.995);
        %slinind=slinind(mask);
        %linind=linind(mask);
        ind=smallu==0;
        energy=energy+sum(newu(ind))/sqrt(dt);

end



end