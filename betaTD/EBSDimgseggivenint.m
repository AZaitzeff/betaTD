function [mapall,dict,energy,gsizes]=EBSDimgseggivenint(mapall,dict,kappa,EBSD,CI,beta,fid,dt,scale)
    numsub=400;
    dx=1/100;
    dy=dx*scale;
    dtstop=2^(-10);
    between=2;
    nt=8;
    flag=0;
    energy=inf;
    gsizes=0;
    while ~flag
        [mapall,dict,kappa,energy,gsizes,flag]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between,numsub);
        dt=dt/2;
    end
    [newmapall,newdict,newkappa,~]=regionmerger(mapall,EBSD,CI,beta,dict,kappa,fid,dtstop,dx,dy,nt,numsub,5);
    [newmapall,newdict,~,newenergy,newgsizes,newflag]=EBSDMStdfast(newmapall,EBSD,CI,beta,newdict,newkappa,fid,dt,dx,dy,dtstop,nt,between,numsub);
    if newflag && newenergy<energy
        mapall=newmapall;
        dict=newdict;
        energy=newenergy;
        gsizes=newgsizes;

    end
   

end