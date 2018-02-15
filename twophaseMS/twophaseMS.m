function twophaseMS(EBSD,CI,map,fid,section)
    [clusterlist,~,labels] = unique(map);
    current=clusterlist;
    M=numel(current);
    new=[];
    z=1;

    tic;updateregion1;toc;

    [clusterlist,~,labels] = unique(map);
    N=numel(clusterlist);
    current=clusterlist;
    M=numel(current);
    tic;updateregion2;toc;

    for iter=1:5
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
    mapedge=zeros(size(map));
    borders = srm_getborders(map);
    mapedge(borders) = 1;
    save(['mappart' num2str(section) num2str(fid) '.mat'],'map','mapedge');
end