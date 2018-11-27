
%addpath('../../../MATLAB/mtex-5.1.1/')
%startup

addpath('../anglelib/')
load('../data/AFbigEBSD')
EBSD=EBSD(1:3:end,1:3:end,:);
CI=CI(1:3:end,1:3:end);
%%
%load('../data/AFoneEBSD')
%EBSD=EBSD(1:2:end,1:2:end,:);
%CI=CI(1:2:end,1:2:end);

%Ks=[2,2];
%[mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks,200,16);
%[mapall,newdict,newkappa]=regionmerging(mapall,newdict,newkappa,EBSD,CI,Ks,.015,200);

%%

%lams=[300,350];
%lams=[300,400];
lams=[100,200,300];
%dts=[.0005,.001,.0025];
dt=2^-6;
Ks=[10,6];
%Ks=[6,6];
dx=1/100*3;
dy=1/100*3;
enes=[0,1];%[1,2];
%for i=1:2
%    for j=1:4

%load('../data/simEBSD')

%EBSD=EBSD((i-1)*600+1:i*600,(j-1)*600+1:j*600,:);
%CI=CI((i-1)*600+1:i*600,(j-1)*600+1:j*600);

%[alphacolors,oM]=crystalcolormaps(EBSD,1);
[m,n,z]=size(EBSD);
for fid=lams
        for enec=enes
    for ztest =1:4
        
        ztest
tic;
[mapall,dict,kappa]=initializeEBSDfast(EBSD,CI,Ks,20,50,50);
[mapall,dict,kappa]=regionmerging(mapall,dict,kappa,EBSD,CI,Ks,.5);

[mapall,dict,kappa]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,dt,dx,dy,(2^-12),enec);
%[mapall,dict,kappa]=EBSDMStd(mapall,sEBSD,sCI,dict,kappa,fid,2^-8,dx,dy,(2^-12));
toc;


energy=EBSDtdE(mapall,EBSD,CI,dict,fid,(2^-12),dx,dy,enec);

%save(['results/AF' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' num2str(ztest) ],'mapall','dict','kappa','energy')


save(['results/AFbig' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt'  num2str(ztest) ],'mapall','dict','kappa','energy')
%addpath('../../../MATLAB/mtex-5.1.1/')
%startup

% mapedgeall=zeros(size(mapall));
% borders = srm_getborders(mapall);
% mapedgeall(borders) = 1;
% imagesc(alphacolors);hold on
% imagesc(repmat(1-mapedgeall,[1,1,3]),'AlphaData',mapedgeall)
% title(num2str(energy))
% axis equal
% print(['AF' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt' num2str(i) num2str(j) 'part' num2str(ztest), ],'-dpng')
% close
    end
    end
    end
 %   end
%end