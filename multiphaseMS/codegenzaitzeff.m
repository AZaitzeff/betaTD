% %%
% 
% vectorType = coder.typeof(1==1, [2000 2000], [true true]);
% codegen makerowcolmapsz -args {vectorType,(7),(7),(7),(7)}
% 
% %%
% 
% vectorType1 = coder.typeof(1, [2000 2000], [true true]);
% vectorType2 = coder.typeof(1, [4000 801], [false false]);
% vectorType3 = coder.typeof(1, [4000 1], [false false]);
% codegen ADIz -args {vectorType1,1,1,1,vectorType2,vectorType2,vectorType3,vectorType3,7,7}

%%
%vectorType = coder.typeof(1, [1 500], [false true]);
%codegen simthompson.m -args {1,1,vectorType,1}

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
vectorType3 = coder.typeof(1, [K 4], [true false]);
vectorType4 = coder.typeof(1, [K 1], [true false]);
codegen initializeEBSDfast -args {vectorType2,vectorType1,[1 1],1,1,1,1}
codegen regionmerging -args {vectorType1,vectorType3,vectorType4,[1 1],1}
vectorType1f = coder.typeof(1, [M N], [true true]);
vectorType2f = coder.typeof(1, [M N 3], [true true false]);
codegen EBSDMStdfast -args {vectorType1f,vectorType2f,vectorType1f,vectorType3,vectorType4,1,1,1,1,1,1,1}


%[mapall,dict,kappa]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,dt,dx,dy,(2^-10),enec);
%[mapall,dict,kappa]=EBSDMStd(mapall,sEBSD,sCI,dict,kappa,fid,2^-8,dx,dy,(2^-12));