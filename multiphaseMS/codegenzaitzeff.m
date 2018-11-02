%%

vectorType = coder.typeof(1==1, [1000 1000], [true true]);
codegen makerowcolmapsz -args {vectorType,(7),(7),(7),(7)}

%%

vectorType1 = coder.typeof(1, [1000 1000], [true true]);
vectorType2 = coder.typeof(1, [2000 401], [false false]);
vectorType3 = coder.typeof(1, [2000 1], [false false]);
codegen ADIz -args {vectorType1,1,1,1,vectorType2,vectorType2,vectorType3,vectorType3,7,7}

%%
%vectorType = coder.typeof(1, [1 500], [false true]);
%codegen simthompson.m -args {1,1,vectorType,1}