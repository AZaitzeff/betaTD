function [mapall,dict,energy,gsizes,flag]=EBSDimgseggivenint(mapall,dict,kappa,EBSD,CI,beta,fid,dt,dx,dy,numsub)
    dtstop=2^(-10);
    between=2;
    nt=6;
    [mapall,dict,kappa,energy,gsizes,flag]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between,numsub);
    if flag
	    [newmapall,newdict,newkappa,newK]=regionmerger(mapall,EBSD,CI,beta,dict,kappa,fid,dtstop*2,dx,dy,nt,numsub,5);
	    [newmapall,newdict,~,newenergy,newgsizes,newflag]=EBSDMStdfast(newmapall,EBSD,CI,beta,newdict,newkappa,fid,dt,dx,dy,dtstop,nt,between,numsub);
        if newenergy<energy
            mapall=newmapall;
            dict=newdict;
            energy=newenergy;
            gsizes=newgsizes;
            flag=newflag;

        end
    end

end