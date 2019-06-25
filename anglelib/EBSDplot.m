function EBSDplot(filename)
close all
EBSDtemp=load(['../data/' filename 'EBSD.mat']);
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[m,n]=size(CI);
M=round(n/50);
N=round(m/50);
[alphacolors,oM]=crystalcolormaps(EBSD,1);
[betacolors,oM]=crystalcolormaps(EBSD,0);
imagesc(alphacolors);
title('Alpha Phase Orientations','FontSize',16);
f=gca;
f.PlotBoxAspectRatio=[M,N,1];
figure
imagesc(betacolors);
title('Beta Phase Orientations','FontSize',16);
f=gca;
f.PlotBoxAspectRatio=[M,N,1];
figure
imagesc(CI);colormap gray;colorbar
title('Confidence Index','FontSize',16)
f=gca;
f.PlotBoxAspectRatio=[M,N,1];
figure
imagesc(EBSDtemp.betas);colormap gray;colorbar
title('Betas','FontSize',16)
f=gca;
f.PlotBoxAspectRatio=[M,N,1];
end