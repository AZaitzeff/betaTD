function [map,dict,K]=zpilchank(EBSD,IQ,CI,beta,alphatol,mincolsize,vartol,betatol,betatol2,numwanted)
    [m,n,~]=size(EBSD);
    [map,K]=floodfill(EBSD,beta,IQ,CI,alphatol,mincolsize);
    [neigharray,sizeneigh]=  findneigh(map,K);
    total=m*n;
    EBSDf=reshape(EBSD, [total 3]);
    quat=E313toq(EBSDf');
    [dict,~,IQmean,betamean]=estimatemeanfast(quat,CI,IQ,beta,map,K,400);
    [betadict,filled]=Pilchank(neigharray,sizeneigh,dict,betamean,K,vartol,betatol);
    [replacemap,newbetadict,K]=postpilchank(neigharray,sizeneigh,betadict,IQmean,filled,K,betatol2);
    newmap=changemapleavezeros(map,replacemap);
    [map,dict,K]=limitnummap(newmap,newbetadict,K,numwanted);
    
