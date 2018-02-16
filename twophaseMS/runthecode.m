


parpool(4)
parfor section=1:8
    
    EBSDtemp=load('../data/scan4subgrid.mat');
    addpath('../anglelib/')
    mapall=load('initial.mat');
    fid=150;
    factor=50;
    if section==9
        rows=1:300
        cols=1:700
    else
        i=mod(section-1,2)+1;
        j=ceil(section/2);
        rows=(max((i-1)*150-factor+1,1):min((i)*150+factor,300));
        cols=(max((j-1)*175-factor+1,1):min((j)*175+factor,700));
        
    end
    EBSD=EBSDtemp.EBSD(rows,cols,:);
    CI=EBSDtemp.CI(rows,cols);
    map=mapall.map(rows,cols);
    tic;
    twophaseMS(EBSD,CI,map,fid,section);
    toc;

    
end