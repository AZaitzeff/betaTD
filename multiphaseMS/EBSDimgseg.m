function [mapall,dict,energy,gsizes]=EBSDimgseg(EBSD,CI,beta,fid,dt,dx,dy,nr,nc,part)
    dtstop=2^(-10);
    between=2;
    nt=6;
    [mapall,dict,kappa,mid]=initializeEBSDfast(EBSD,CI,beta,nr,nc);
    [mapall,dict,kappa]=regionmerging(mapall,dict,kappa,nr,nc,mid,fid);
    if part
        K=size(dict,1);
        if K<200
            [mapall,dict,~,energy,gsizes]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between);
        else
            energy=inf;
            gsizes=400;
        end
    else
        [mapall,dict,~,energy,gsizes]=EBSDMStdfast(mapall,EBSD,CI,beta,dict,kappa,fid,dt,dx,dy,dtstop,nt,between);
    end
end