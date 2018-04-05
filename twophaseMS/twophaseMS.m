function twophaseMS(EBSD,CI,betas,map,fid,betathres,section,filename)
    subnum=500;
    betas(CI<betathres)=0;
    dict=estimatebetas(EBSD,CI,betas,map,[],0,1,subnum);
    [clusterlist,~,labels] = unique(map);
    current=clusterlist;
    M=numel(current);
    new=[];
    z=1;

    updateregion1;
    if numel(new)>0
        dict=estimatebetas(EBSD,CI,betas,map,new,dict,1,subnum);
    end
    [clusterlist,~,labels] = unique(map);
    N=numel(clusterlist);
    current=clusterlist;
    M=numel(current);
    updateregion2;

    for iter=1:5
         [current,~,labels] = unique(new);
         M=numel(current);
         new=[];
         z=1;
         updateregion1
         if numel(new)>0
            dict=estimatebetas(EBSD,CI,betas,map,new,dict,1,subnum);
         end
         [clusterlist,~,labels] = unique(map);
         N=numel(clusterlist);
         [current,~,labels] = unique(new);
         M=numel(current);
        updateregion2
        
    end
    cleanup;
    save(['results/' filename 'part' num2str(section) num2str(fid) '.mat'],'map','dict');
end