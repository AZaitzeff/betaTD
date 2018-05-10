filename='test1';
%run2pMS(filename,'beta0',[50,100,150,200,250,300],500,4,0);
%run2pMS(filename,'beta25',[50,100,150,200,250,300],500,4,.25);
%run2pMS(filename,'beta50',[50,100,150,200,250,300],500,4,.5);
%run2pMS(filename,'beta75',[50,100,150,200,250,300],500,4,.75);
run2pMS(filename,'beta1',[250,500,750],500,4,1);


filename='test2';
%run2pMS(filename,'betat250',[50,100,150,200],400,4,.5);
run2pMS(filename,'betat21',[200],400,4,1);

%filename='RX';
%run2pMS(filename,'beta50RX',[100,150,200],400,4,.5);
%run2pMS(filename,'betat1RX',[100,150,200],400,4,1);

filename='sim';
run2pMS(filename,'betatsim1',[50,100,150,200],300,4,1);