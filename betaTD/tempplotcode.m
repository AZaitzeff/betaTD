addpath('../anglelib/')
%addpath('../../../MATLAB/mtex-5.1.1/')
%startup
%b                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               addpath('../anglelib/')
EBSDtemp=load(['../data/' 'RX' 'EBSD.mat']);
EBSD=EBSDtemp.EBSD;
beta=logical(EBSDtemp.betas);
[M,N]=size(EBSDtemp.CI);
load('results/RXpatchfull.mat')
%load('../floodfill/results/RXff.mat')
% [colorsbnd,bnds]=  bndmap(mapall,dict,[]);
% betaEBSD=converttobetamap(EBSD,beta,dict,mapall);
%[colors,oM]=crystalcolormaps(betaEBSD,0);
[colors,oM]=crystalcolormaps(truebetaEBSD,0);
%[colors,oM]=crystalcolormaps(EBSD,1);
imagesc(colors);hold on
%imagesc(colorsbnd,'AlphaData',bnds)

f=gca;
f.PlotBoxAspectRatio=[ceil(N/50),ceil(M/50),1];

print(['figures/' 'RX' 'patch'],'-dpng')
close