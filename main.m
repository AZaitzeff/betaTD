addpath('pilchank/')
addpath('betaTD/')
addpath('anglelib/')

% sizepar should be:
% 1: beta grains around 50 by 50 pixels or smaller
% 2: beta grains around 100 by 100 pixels (a good pick if you do not know)
% 3: beta grains around 200 by 200 pixels or larger
sizepar=2; %change
mex=0; % If have matlab coder, you can change this to 1.

runbetatdpilchank(nickname,nickname,sizepar,mex);

%runbetatd(nickname,nickname,2,10,2^-4,100,30); %Different method that uses
%random initialization. 

%plots the results

set(0,'DefaultFigureVisible','off');
thres=[];%highlights bnds under certain tolorances (degrees) e.g. [2,5,10]


if sizepar==1
    fids=[100,150,200,250];
elseif sizepar==2
    fids=[75,100,125,150,175];
else
    fids=[25,50,75,100];
end

for fid=fids
    plottingcode(nickname,[nickname num2str(fid)],thres,1,1);
end

plottingcode(nickname,[nickname num2str(fid)],[],1,3);


set(0,'DefaultFigureVisible','on');