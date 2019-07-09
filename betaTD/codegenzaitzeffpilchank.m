function codegenzaitzeffpilchank(M,N)
addpath('../othermethods')
vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);
%vectorType3 = coder.typeof(1, [K 4], [true false]);
%vectorType4 = coder.typeof(1, [K 1], [true false]);
%codegen regionmerging -args {vectorType1,vectorType3,vectorType4,[1 1],1}
%codegen EBSDMStdfast -args {vectorType1,vectorType2,vectorType1,vectorType1l,vectorType3,vectorType4,1,1,1,1,1,1,1}
codegen zpilchank -args {vectorType2,vectorType1,vectorType1,vectorType1l,1,1,1,1,1}
end
