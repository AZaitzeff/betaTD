load('../data/scan4subgrid.mat')
addpath('../anglelib/')
load('initial.mat')
EBSD=EBSD(31:130,71:170,:);
CI=CI(31:130,71:170);
map1=map(31:130,71:170);
for fid=[75,100,125,150,200]
    %EBSD=EBSD(1:300,1:700,:);
    %CI=CI(1:300,1:700);
    map=map1;
    


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
        save(['mapf' num2str(iter) num2str(fid) '.mat'],'map');
         [clusterlist,~,labels] = unique(map);
         N=numel(clusterlist);
         [current,~,labels] = unique(new);
         M=numel(current);
        updateregion2
        save(['map' num2str(iter) num2str(fid) '.mat'],'map');
    end
end