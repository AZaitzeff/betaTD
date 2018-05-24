filename='test1';
%run2pMS(filename,'beta0',[50,100,150,200,250,300],500,4,0);
%run2pMS(filename,'beta25',[50,100,150,200,250,300],500,4,.25);
%run2pMS(filename,'beta50',[50,100,150,200,250,300],500,4,.5);
%run2pMS(filename,'beta75',[50,100,150,200,250,300],500,4,.75);
run2pMS(filename,'beta1',[200,400,800,1600],300,4,1);


filename='test2';
%run2pMS(filename,'betat250',[50,100,150,200],400,4,.5);
run2pMS(filename,'beta21',[200,400,800,1600],300,4,1);

%filename='RX';
%run2pMS(filename,'beta50RX',[100,150,200],400,4,.5);
%run2pMS(filename,'beta1RX',[400,600,800],400,4,1);

%filename='sim';
%run2pMS(filename,'betatsim1',[200,300,400,500],300,4,1);