function plotebsd(EBSD,valued,beta)


[m,n,~]=size(EBSD);
total=m*n;
valued=valued(:);
alphacolors=zeros(total,3);
EBSDflat=reshape(EBSD,[total,3]);
[partcolors,oM]=crystalcolormaps(EBSDflat(valued,:),~beta);
alphacolors(valued,:)=partcolors;
alphacolors=reshape(alphacolors, [m,n,3]);
M=round(n/50);
N=round(m/50);

imagesc(alphacolors);
%title('Alpha Phase Orientations','FontSize',16);
f=gca;
f.PlotBoxAspectRatio=[M,N,1];
end

