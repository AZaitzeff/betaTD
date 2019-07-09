addpath('pilchank/')
addpath('betaTD/')

%change
fids=[25,50,100];
runMStdpilchank(nickname,nickname,2^-4,fids,0,200);


%plots the results
set(0,'DefaultFigureVisible','off');

