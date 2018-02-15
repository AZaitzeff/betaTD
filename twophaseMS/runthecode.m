
for section=1:8
    EBSDtemp=load('../data/scan4subgrid.mat');
    addpath('../anglelib/')
    mapall=load('initial.mat');
    fid=150;
    factor=15;
    i=mod(section-1,2)+1;
    j=ceil(section/2);
    rows=(max(i*150-factor,0):min((i+1)*150+factor,300));
    cols=(max(j*175-factor,0):min((j+1)*175+factor,700));
    EBSD=EBSDtemp.EBSD(rows,cols,:);
    CI=EBSDtemp.CI(rows,cols);
    map=mapall.map(rows,cols);
    


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