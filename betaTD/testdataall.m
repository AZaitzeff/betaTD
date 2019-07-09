filesave='RX';
fids=[100,150,200,250,300,350,400];

numfids=numel(fids);
score=zeros(1,numfids);
score2=zeros(1,numfids);
gsizes=zeros(1,numfids);
runcheck=20;
load(['../data/' filesave 'EBSD.mat'])
beta=logical(betas);
for z=1:numfids
    fid=fids(z);
    [~,I]=min(energies);
    var=load(['results/' filesave num2str(round(fid))]);
    vals=matchmetric(var.mapall,var.dict);
    [val,~]=fidmetric(EBSD,CI,beta,var.mapall,var.dict);
    score2(z)=val;
    %figure
    %hist(vals)
    if isempty(vals)
         score(z)=60;
    else
        score(z)=prctile(vals,1);
    end
    %mapall=var.mapall;
    %save(['results/' filesave 'iter' num2str(round(fid))],'mapall');
end