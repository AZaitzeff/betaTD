load('../data/scan4subgrid.mat')
addpath('../anglelib/')
load('initial.mat')
for fid=[100,150,200,250]
    %EBSD=EBSD(1:300,1:700,:);
    %CI=CI(1:300,1:700);
    EBSD=EBSD(31:180,31:180,:);
    CI=CI(31:180,31:180);
    map=map(31:180,31:180);


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