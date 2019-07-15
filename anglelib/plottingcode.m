function plottingcode(filename,filesave,thres,step,which)
if nargin<4
    step=1;
end
if nargin<5
    which=1;
end

if ~exist(['figures/' filename],'dir')
    status = mkdir(['figures/' filename]);
end


EBSDtemp=load(['data/' filename 'EBSD.mat']);
EBSD=EBSDtemp.EBSD(1:step:end,1:step:end,:);
CI=EBSDtemp.CI(1:step:end,1:step:end);
[m,n]=size(CI);
M=round(n/50);
N=round(m/50);
total=m*n;
alphacolors=zeros(total,3);
beta=reshape(logical(EBSDtemp.betas(1:step:end,1:step:end)), [total,1]);
EBSDflat=reshape(EBSD,[total,3]);
if which==0 || which==2 || which==3
    if sum(~beta)>0
        [alphacolorstemp,oM]=crystalcolormaps(EBSDflat(~beta,:),1);
        alphacolors(~beta,:)=alphacolorstemp;
    end
    if sum(beta)>0
        [betacolorstemp,oM]=crystalcolormaps(EBSDflat(beta,:),0);
        alphacolors(beta,:)=betacolorstemp;
    end
    
    

    alphacolors=reshape(alphacolors, [m,n,3]);
end


if which==0 || which==3
    imagesc(alphacolors);
    %title('Alpha Phase Orientations','FontSize',16);
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    print(['figures/' filename '/' filename 'Alphadata'],'-dpng')
    close
    imagesc(CI);colormap gray;colorbar
    %title('Confidence Index','FontSize',16)
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    print(['figures/' filename '/' filename 'CI'],'-dpng')
    close
    imagesc(EBSDtemp.betas);colormap gray;colorbar
    %title('Confidence Index','FontSize',16)
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    print(['figures/' filename '/' filename 'betas'],'-dpng')
    close
end

Z=numel(thres);
colormapz=[[0,0,0];[0.9290, 0.6940, 0.1250];[0.8500, 0.3250, 0.0980];...
    [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];

clear EBSD CI EBSDtemp;

maps=load(['results/' filesave]);
mapall=maps.mapall;
betaEBSD=maps.betaEBSD;
energy=maps.energy;
dict=maps.dict;
[colorsbnd,bnds]=  bndmap(mapall,dict,thres);
if which==0 || which==1
    [colors,oM]=crystalcolormaps(betaEBSD,0);
    imagesc(colors);hold on
    imagesc(colorsbnd,'AlphaData',bnds)
    %title(['Beta Phase Orientations, Fidelity ' num2str(fid) ', Energy ' num2str(round(energy))],'FontSize',12);
    if Z>0
        h = zeros(Z, 1);
        hname=cell(Z, 1);
        for i=1:Z
            h(i) = plot(NaN,NaN,'Color',colormapz(i+1,:));
            hname{i}=num2str(round(thres(i),1));

        end
        legend(h, hname,'Location','NorthEastOutside');
    end
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    print(['figures/' filename '/' filesave 'Betadata'],'-dpng')
    close
end

if which==0 || which==2
    imagesc(alphacolors);hold on
    imagesc(colorsbnd,'AlphaData',bnds)
    %title(['Beta Grains Overlayed on Alpha Orientation Data Fidelity ' num2str(fid)],'FontSize',12);
    f=gca;
    f.PlotBoxAspectRatio=[M,N,1];
    print(['figures/' filename '/' filesave 'Betaoverlay'],'-dpng')
    close
end


