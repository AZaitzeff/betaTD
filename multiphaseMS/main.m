%names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
% for i=5:8
%     filename=names{i};
%     name=names{i};
%     runMStd(filename,name,12,24,50);
% end
% 
% 
%  names={'rand','sim2f','sim3f','sim4f','AFbig'};
%  for i=1:5
%      filename=names{i};
%      name=names{i};
%      runMStd(filename,name,12,24,50);
%  end


sizes=[200,100,50,25,12.5,6.25];
 for gs=sizes
     filename='AFone';
     name=['AFone' num2str(gs) 's'];
     runMStd(filename,name,12,12,gs);
 end