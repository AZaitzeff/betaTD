function twophaseMS(EBSD,CI,map,fid,part,filename)
    [clusterlist,~,labels] = unique(map);
    current=clusterlist;
    M=numel(current);
    new=[];
    z=1;

    updateregion1;

    [clusterlist,~,labels] = unique(map);
    N=numel(clusterlist);
    current=clusterlist;
    M=numel(current);
    updateregion2;

    for iter=1:5
        iter
         [current,~,labels] = unique(new);
         M=numel(current);
         new=[];
         z=1;
         updateregion1
         [clusterlist,~,labels] = unique(map);
         N=numel(clusterlist);
         [current,~,labels] = unique(new);
         M=numel(current);
        updateregion2
        
    end
    cleanup;
    save(['results/' filename 'part' num2str(part) num2str(fid) '.mat'],'map');
end