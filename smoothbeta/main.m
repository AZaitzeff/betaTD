load('testEthan.mat');
addpath('../anglelib/');
[M,N,~]=size(EBSD);

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);

codegen smoothonegrain -args {vectorType2,vectorType1,vectorType1l,1,1,1,vectorType1}


%smooth=smoothonegrain_mex(EBSD,CI,beta,nt,dt,fid,w);