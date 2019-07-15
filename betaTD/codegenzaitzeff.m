function codegenzaitzeff(M,N)

vectorType1 = coder.typeof(1, [M N], [false false]);
vectorType1l = coder.typeof(1==1, [M N], [false false]);
vectorType2 = coder.typeof(1, [M N 3], [false false false]);

codegen EBSDimgseg -args {vectorType2,vectorType1,vectorType1l,1,1,1,1,1}
end
