 names={'sim','AFone','RX','rand','mapcenter','mapedge','hardbot','hardmid','hardtop','sim2','AFbeta','AFbig'};
fids=100:50:300;
numfids=numel(fids);
 score=zeros(1,numfids);
for i=[5]
name=names{i};
for z=1:numfids

    vars=load(['results/' name num2str(round(z*50+50))]);
    [val]=matchmetric(vars.mapall,vars.dict);
    score(z)=prctile(val,1);
end

plot(fids,score) 
print(['figures/' name '/' name 'match'],'-dpng');
close

I=find(score<1,1)-1;
if isempty(I)
    fid=fids(numfids);
elseif I==0
    fid=100;

else
    fid=fids(I);
end
 
%%

EBSDtemp=load(['../data/' name 'EBSD.mat']);
EBSD=EBSDtemp.EBSD;
CI=EBSDtemp.CI;
[M,N]=size(CI);
num=20;
beta=logical(EBSDtemp.betas);
fidterms=zeros(1,num);
for n=1:num
    vars=load(['results/' name num2str(round(fid)) num2str(n)]);
    [val]=fidmetric(EBSD,CI,beta,vars.mapall,vars.dict,1);
    fidterms(n)=val;
end

[I,score,probmap]=probmetric(['results/' name num2str(round(fid))],5,num);
[I,conval,conmap]=confidencemap2(['results/' name num2str(round(fid))],M,N,num);
name
vars=load(['results/' name num2str(round(fid))]);
1/vars.conval
conval
1/var(fidterms)
score
end