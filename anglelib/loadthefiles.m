addpath('../../../MATLAB/mtex-5.1.1/')
startup

Files=dir('../../Ti64/');
filenames={'AFone','AFbig','AFalpha','RX','hardbot','hardmid','hardtop','mapcenter','mapedge'};
%names={'2B-B-1-1HT1-RF-cleaned.ang',};
for i=[3:4 6:11]
    [EBSD,CI,IQ,betas]=ang2mat(fileload,filenames{i-2});
end