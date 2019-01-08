names={'sim','AFone','RX','rand','AFbig','mapcenter','mapedge','hardbot','hardmid','hardtop'};
for i=2:2
   for fid=25:25:300
       filename=names{i};
       name=names{i};
       runMStdsimple(filename,name,16,16,2^-5,fid,50);
   end
end
 
% names={'sim','AFone','RX','rand','AFbeta','mapcenter','mapedge','hardbot','hardmid','hardtop','AFbig'};
% for i=1:11
%   filename=names{i};
%   name=names{i};
%   runMStd(filename,name,16,32,1);
% end
% 
%names={'sim','sim2f','sim3f','sim4f'};
%for i=1:4
%  filename=names{i};
%  name=names{i};
%  runMStdsimple(filename,name,12,24,2^-5,100,50);
%end

%  names={'sim','sim2f','sim3f','sim4f'};
%  names2={'sim1f','sim2f','sim3f','sim4f'};
%  for i=2:4
%      filename=names{i};
%      name=names2{i};
%      runMStd(filename,name,1,10,50);
%  end

%  for KSfac=[1,2,3,4]
%      filename='AFone';
%      name=['AFone' num2str(KSfac) 'fac'];
%      runMStd(filename,name,12,24,50,KSfac);
%  end
 