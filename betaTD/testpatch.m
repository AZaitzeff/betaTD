names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
addpath('../anglelib/')
ind=2;
name=names{ind}; 
load(['results/' name 'scores']);

fids=[50,100,150,200,250,300];
num=9;
test=zeros(6,10);
for z=1:6
    for some=1:10
    ind=randsample(9,num,'false');
    threeene=energies(ind,z);
    scores=score3(ind,z);
    [~,I]=min(threeene);
    thescore=scores(I);
    test(z,some)=thescore;
    end
    
    
end
plot(fids,test);
