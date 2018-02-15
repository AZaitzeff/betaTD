parpool(5)
parfor section=1:8
    EBSDtemp=load('../data/scan4subgrid.mat');
    addpath('../anglelib/')
    mapall=load('initial.mat');
    fid=150;
    factor=15;
    i=mod(section-1,2)+1;
    j=ceil(section/2);
    rows=(max(i*150-factor,0):min((i+1)*150+factor,300));
    cols=(max(j*175-factor,0):min((j+1)*175+factor,700));
    EBSD=EBSDtemp.EBSD(rows,cols,:);
    CI=EBSDtemp.CI(rows,cols);
    map=mapall.map(rows,cols);
    twophaseMS(EBSD,CI,map,fid,section);


    
end