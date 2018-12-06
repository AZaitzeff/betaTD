function codegenzaitzeff(M,N,K)

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
vectorType3 = coder.typeof(1, [K 4], [true false]);
vectorType4 = coder.typeof(1, [K 1], [true false]);
codegen initializeEBSDfast -args {vectorType2,vectorType1,vectorType1l,[1 1],1,1,1,1}
codegen regionmerging -args {vectorType1,vectorType3,vectorType4,[1 1],1}
codegen EBSDMStdfast -args {vectorType1,vectorType2,vectorType1,vectorType1l,vectorType3,vectorType4,1,1,1,1,1,1,1}

end
%[mapall,dict,kappa]=EBSDMStd(mapall,EBSD,CI,dict,kappa,fid,dt,dx,dy,(2^-10),enec);
%[mapall,dict,kappa]=EBSDMStd(mapall,sEBSD,sCI,dict,kappa,fid,2^-8,dx,dy,(2^-12));