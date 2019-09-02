addpath('pilchak/')
addpath('betaTD/')
addpath('anglelib/')

% sizepar should be:
% 1: beta grains around 50 by 50 pixels or smaller
% 2: beta grains around 100 by 100 pixels (a good pick if you do not know)
% 3: beta grains around 200 by 200 pixels or larger
sizepar=2; %change
mex=0; % If have matlab coder, you can change this to 1.

%nickname='AF15';
%runbetatdpilchak(nickname,nickname,sizepar,mex);


fids=[75,100,125,150];
for fid=fids
    runbetatd(nickname,nickname,4,20,fid,2^-4,30); %Different method that uses
end
%random initialization. 


