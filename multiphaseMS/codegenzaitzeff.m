%%

vectorType = coder.typeof(1==1, [2000 2000], [true true]);
codegen makerowcolmapsz -args {vectorType,(7),(7),(7),(7)}

%%

vectorType1 = coder.typeof(1, [2000 2000], [true true]);
vectorType2 = coder.typeof(1, [4000 801], [false false]);
vectorType3 = coder.typeof(1, [4000 1], [false false]);
codegen ADIz -args {vectorType1,1,1,1,vectorType2,vectorType2,vectorType3,vectorType3,7,7}

%%
%vectorType = coder.typeof(1, [1 500], [false true]);
%codegen simthompson.m -args {1,1,vectorType,1}

vectorType1 = coder.typeof(1, [400 400], [false false]);
vectorType2 = coder.typeof(1, [400 400 3], [false false false]);
vectorType3 = coder.typeof(1, [100 4], [true false]);
vectorType4 = coder.typeof(1, [100 1], [true false]);
codegen EBSDMStdfast -args {vectorType1,vectorType2,vectorType1,vectorType3,vectorType4,1,1,1,1,1,1,1}