load('../data/scan4subgrid.mat')
addpath('../anglelib/')
load('initial.mat')
fid=50;
EBSD=EBSD(1:300,1:700,:);
CI=CI(1:300,1:700);


[clusterlist,~,labels] = unique(map);
current=clusterlist;
M=numel(current);
new=[];
z=1;

tic;updateregion1;toc;
e=cputime-t

[clusterlist,~,labels] = unique(map);
N=numel(clusterlist);
current=clusterlist;
M=numel(current);
tic;updateregion2;toc;
e=cputime-t

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