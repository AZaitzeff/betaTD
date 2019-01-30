
%addpath('../../../MATLAB/mtex-5.1.1/')
%startup

 addpath('../anglelib/')
 load('../data/AFoneEBSD')
% EBSD=EBSD(1:3:end,1:3:end,:);
% CI=CI(1:3:end,1:3:end);
%%
load('../data/AFbigEBSD')
[M,N]=size(CI);
midm=ceil(M/2);
midn=ceil(N/2);
lm=max(midm-200,1);
um=min(midm+200,M);
ln=max(midn-200,1);
un=min(midn+200,N);
suby=lm:um;
subx=ln:un;
M=um-lm+1;
N=un-ln+1;
EBSD=EBSD(suby,subx,:);
CI=CI(suby,subx);
beta=logical(betas(suby,subx));
% EBSD=EBSD(401:800,401:800,:);
% CI=CI(401:800,401:800);
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
nr=10;
nc=10;
codegenzaitzeff(4000,4000);
fid=200;
dt=2^-5;
%% 
for ztest =1:2
        



tic;
[mapall,dict,energy,gsizes]=EBSDimgseg_mex(EBSD,CI,beta,fid,dt,dx,dy,nr,nc);
toc;
% [vals]=matchmetric(mapall,dict);
% prctile(vals,1)
%     %mapall=var.mapall;
% round(prctile(gsizes,5))

%save(['results/AFbeta' num2str(fid) 'fid' num2str(ztest) ],'mapall','dict','gsizes','energy')

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