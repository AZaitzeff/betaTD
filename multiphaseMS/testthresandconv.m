
%addpath('../../../MATLAB/mtex-5.1.1/')
%startup

 addpath('../anglelib/')
% load('../data/AFbigEBSD')
% EBSD=EBSD(1:3:end,1:3:end,:);
% CI=CI(1:3:end,1:3:end);
%%
load('../data/AFoneEBSD')
%EBSD=EBSD(1:2:end,1:2:end,:);
%CI=CI(1:2:end,1:2:end);
%load('../data/simEBSD')
%Ks=[2,2];
%[mapall,newdict,newkappa]=initializeEBSDfast(EBSD,CI,Ks,200,16);
%[mapall,newdict,newkappa]=regionmerging(mapall,newdict,newkappa,EBSD,CI,Ks,.015,200);

%%

%lams=[300,350];

dx=1/100;
dy=1/100;
%[1,2];
%for i=1:2
%    for j=1:4

%load('../data/simEBSD')

%EBSD=EBSD((i-1)*600+1:i*600,(j-1)*600+1:j*600,:);
%CI=CI((i-1)*600+1:i*600,(j-1)*600+1:j*600);

%[alphacolors,oM]=crystalcolormaps(EBSD,1);
[M,N,z]=size(EBSD);
%nr=10;
%nc=10;
nr=20;
nc=20;
codegenzaitzeff(M,N);
fid=150;
dt=2^-5;
for ztest =1:1
        



tic;
[mapall,dict,energy,gsizes]=EBSDimgseg_mex(EBSD,CI,logical(betas),fid,dt,dx,dy,nr,nc);
toc;


save(['results/AF' num2str(fid) 'fid' num2str(ztest) ],'mapall','dict','gsizes','energy')

%save(['results/sim' num2str(fid) 'fid' num2str(ztest) ],'mapall','dict','kappa','energy')
%save(['results/AFbig' num2str(enec) 'e' num2str(fid) 'fid' num2str(round(dt*100)) 'dt'  num2str(ztest) ],'mapall','dict','kappa','energy')
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
 %   end
%end