
codegen VMFEMfast -args {ones(200,4),ones(4,4,144),1,ones(1,4),1}

vectorType = coder.typeof(1, [200 4], [true false]);

codegen VMFEMzfast -args {vectorType,ones(4,4,144),vectorType,ones(4,4,24),1,ones(1,4),1}