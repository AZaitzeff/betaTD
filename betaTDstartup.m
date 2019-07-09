%CHANGE
addpath('../../MATLAB/mtex-5.1.1/') %put path to mtex distribtuion here
startup;

addpath('anglelib/')


CS = {... 
  'notIndexed',...
  crystalSymmetry('622', [2.95 2.95 4.68], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Titanium (Alpha)', 'color', 'light blue'),...
  crystalSymmetry('432', [3.31 3.31 3.31], 'mineral', 'Titanium (Beta)', 'color', 'light green')};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names (CHANGE)

%user created nickname
nickname= 'scanA'; 

% path to files
pname = '/Users/azaitzeff/Documents/Research/Ti64';

% which file to be imported
fname = [pname '/Ti64_1100C_0p5hr_scan20190525A']; 
%% Import the Data

% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');
[EBSD,CI,IQ,betas,scale]=ang2mat(ebsd,0); %if 'Titanium - Beta' change to 0 to 1

save(['data/' nickname 'EBSD'],'CI','IQ','EBSD','betas','scale')% puts the code in a format
