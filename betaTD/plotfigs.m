names={'AFbig','R1','R2','AFone','RX','rand','sim','AFbeta','AF225','AFpancake','fine'};


addpath('../anglelib/')

addpath('../../../MATLAB/mtex-5.1.1/')
startup

set(0,'DefaultFigureVisible','off');


%plottingcode('AFbigrow2col2','AFbigrow2col250',[],1,1);
%plottingcode('AFbigrow2col2','AFbigrow2col2100',[],1,1);

 %plottingcode('Merged','Merged12',[],1,1);
 plottingcode('Merged','Merged25',[],1,1);
% plottingcode('Merged','Merged40',[],1,1);
% plottingcode('Merged','Merged50',[],1,1);

% 
%filename='sim2';
%name=['sim225'];
%plottingcode(filename,name,[],1,0);

 
% filename='mapcenter';
% name=['mapcenter' 'beta'];
% plottingcode(filename,name,[],1,0);
% for fid=50:50:400
% filename='sim';
% name=['sim' num2str(fid)];
% plottingcode(filename,name,[],1,0);
% end
% i=3;
% filename=names{i};
% name=[names{i} num2str(100)];
% plottingcode(filename,name,[],6,3);

%  z=1;
%  for i=2:2
%     for j=2:2
%     filename=[names{z} 'row' num2str(2) 'col' num2str(2)];
%     name=[filename '8gs5dt200'];
%     plottingcode('AFbig',name,[],1,1);
%     end
%  end
%  z=10;
%  for i=1:4
%     for j=1:4
%              filename=[names{z} 'row' num2str(i) 'col' num2str(j)];
%              name=filename;
%              plottingcode(filename,[name '50'],[],1,1);
%              %runMStdsimple(filename,name,5,20,-1,150,12,200,0);
%     end
%  end


%plottingcode(filename,name(1:end-4),[],1,3);

% end

% for i=[4]
% 
% filename=names{i};
% 
% S = dir(['results/' filename '*.mat']);
% 
% for file=S'
% name=file.name;
% plottingcode(filename,name(1:end-4),[],1,1);
% end
% % for fid=100:50:300
% % filename=names{i};
% % name=[names{i} num2str(fid)];
% % plottingcode(filename,name,[],1,1);
% % end
% 
% end
%     filename='AFbig';
%     name=['AFbig150'];
%     plottingcode(filename,name,[],1,0);

% 
%  for enec=[1,2]
%          if enec>=0
%              name=['AFone' num2str(10*enec)];
%          else
%              name=['AFonen' num2str(-10*enec)];
%          end
%          filename='AFone';
%          plottingcode(filename,name,[100,200],[1,2,5,10],2);
%  end
% % 
%  for enec=[0]
%          if enec>=0
%              name=['AFbig' num2str(10*enec)];
%          else
%              name=['AFbign' num2str(-10*enec)];
%          end
%          filename='AFbig';
%          plottingcode(filename,name,[200],[],4);
%  end

set(0,'DefaultFigureVisible','on');