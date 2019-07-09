addpath('pilchank/')
addpath('betaTD/')

%Good defualts more detail to be added
fids=[25,50,100];
runMStdpilchank(nickname,nickname,2^-4,fids,0,200);


%plots the results

set(0,'DefaultFigureVisible','off');
thres=[];%highlights bnds under certain tolorances (degrees) e.g. [2,5,10]
for fid=fids
    plottingcode(nickname,[nickname num2str(fid)],thres,1,1);
end

plottingcode(nickname,[nickname num2str(fid)],[],1,3);


set(0,'DefaultFigureVisible','on');