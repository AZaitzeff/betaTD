function [map,dict,K]=zpilchak(EBSD,IQ,CI,beta,alphatol,mincolsize,vartol,betatol,betatol2,numwanted)
    [m,n,~]=size(EBSD);
    [map,K]=floodfill(EBSD,beta,IQ,CI,alphatol,mincolsize);
    [neigharray,sizeneigh]=  findneigh(map,K);
    total=m*n;
    EBSDf=reshape(EBSD, [total 3]);
    quat=E313toq(EBSDf');
    [dict,~,IQmean,betamean]=estimatemeanfast(quat,CI,IQ,beta,map,K,400);
    [betadict,filled]=Pilchak(neigharray,sizeneigh,dict,betamean,K,vartol,betatol);
    [replacemap,newbetadict,K]=postpilchak(neigharray,sizeneigh,betadict,IQmean,filled,K,betatol2);
    newmap=changemapleavezeros(map,replacemap);
    [map,dict,K]=limitnummap(newmap,newbetadict,K,numwanted);
    
