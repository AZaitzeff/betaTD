names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
for i=6:8
    filename=names{i};
    name=names{i};
    runMStd(filename,name,12,24,50);
end


 names={'rand','sim2f','sim3f','sim4f','AFbig'};
 for i=1:5
     filename=names{i};
     name=names{i};
     runMStd(filename,name,12,24,50);
 end

%  names={'sim','sim2f','sim3f','sim4f'};
%  names2={'sim1f','sim2f','sim3f','sim4f'};
%  for i=2:4
%      filename=names{i};
%      name=names2{i};
%      runMStd(filename,name,1,10,50);
%  end

% sizes=[200,100,50,25,12.5,6.25];
%  for gs=sizes
%      filename='AFone';
%      name=['AFone' num2str(gs) 's'];
%      runMStd(filename,name,12,12,gs);
%  end