
addpath('../../../MATLAB/mtex-5.1.1/')
startup



addpath('../anglelib/')
%%
%load('../data/AFoneEBSD')
load('../data/simEBSD')
% EBSD=EBSD(1:2:end,1:2:end,:);
% CI=CI(1:2:end,1:2:end);
%Ks=[2,2];
%[mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks,200,16);
%[mapall,newdict,newkappa]=regionmerging(mapall,newdict,newkappa,EBSD,CI,Ks,.015,200);

%%
thres=[];
Z=numel(thres);
colormapz=[[0,0,0];[0.9290, 0.6940, 0.1250];[0.8500, 0.3250, 0.0980];...
    [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];
%lams=[200,300];
%lams=[50,100,150,200,250,300];
lams=[200];
%dts=[.0005,.001,.0025];
dt=2^-6;
Ks=[10,6];
%Ks=[6,6];
dx=1/100;
dy=1/100;
enes=[0];%[1,2];



%  load('../data/randEBSD')
  
   %load('../data/simEBSD.mat')
%  EBSD=EBSD(1:600,1:600,:);
%  CI=CI(1:600,1:600);
%EBSD=EBSD(1:3:end,1:3:end,:);
% CI=CI(1:3:end,1:3:end);
[alphacolors,oM]=crystalcolormaps(EBSD,1);

%%

[m,n,z]=size(EBSD);
M=round(n/50);
N=round(m/50);
for fid=25:25:300
        


%load(['results/AFone' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' num2str(ztest), ])
%load(['results/AF' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' num2str(ztest)])
%load(['results/AF' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' '11part' num2str(ztest)])
%load(['results/AFone' num2str(fid) num2str(ztest) ])
load(['results/mapcenter' num2str(fid) ])
% mapedgeall=zeros(size(mapall));
% borders = srm_getborders(mapall);
% mapedgeall(borders) = 1;
% imagesc(alphacolors);hold on
% imagesc(repmat(1-mapedgeall,[1,1,3]),'AlphaData',mapedgeall)
% title(num2str(energy))
    [colorsbnd,bnds]=  bndmap(mapall,[],[]);
    %[colorsbnd,bnds]=  bndmap(mapall,dict,thres);
    imagesc(alphacolors);hold on
    imagesc(colorsbnd,'AlphaData',bnds)
%     h = zeros(Z, 1);
%     hname=cell(Z, 1);
%     for temp=1:Z
%         h(temp) = plot(NaN,NaN,'Color',colormapz(temp+1,:));
%         hname{temp}=num2str(round(thres(temp),1));
% 
%     end
%     legend(h, hname,'Location','NorthEastOutside');
    %title(num2str(round(energy)))
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    
    
    % axis equal
    print(['mc' num2str(fid)],'-dpng')
    %print(['rand' num2str(fid) 'fid' num2str(ztest), ],'-dpng')
    %print(['AFbigc' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' num2str(ztest), ],'-dpng')
    %print(['sim'  num2str(fid) 'fid' num2str(ztest), ],'-dpng')
 close
end