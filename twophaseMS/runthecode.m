parpool(4)
filename='test1';
run2pMS(filename,[75,100,125,150],300);
filename='test2';
run2pMS(filename,[75,100,125,150],400);

% parpool(5)
% fids=[75,100,125,150,175];
% parfor section=1:8
%     for fid=fids
%         EBSDtemp=load('../data/scan4subgrid.mat');
%         addpath('../anglelib/')
%         mapall=load('initial.mat');
%         factor=25;
%         if section==9
%             rows=1:300
%             cols=1:700
%         else
%             i=mod(section-1,2)+1;
%             j=ceil(section/2);
%             rows=(max((i-1)*150-factor+1,1):min((i)*150+factor,300));
%             cols=(max((j-1)*175-factor+1,1):min((j)*175+factor,700));
% 
%         end
%         EBSD=EBSDtemp.EBSD(rows,cols,:);
%         CI=EBSDtemp.CI(rows,cols);
%         map=mapall.map(rows,cols);
%         tic;
%         twophaseMS(EBSD,CI,map,fid,section);
%         toc;
% 
%     end
% end
% 
% puttogether(fids)