 names={'sim','AFone','RX','rand','mapcenter','mapedge','hardbot','hardmid','hardtop','sim2','AFbeta','AFbig'};
 for i=[12]
    for fid=[25,50,100]
        filename=names{i};
        name=names{i};
        runMStdsimple(filename,name,16,20,2^-4,fid,30);
    end
 end
 
 
% for i=12:12
%     for fid=100:50:400
%         filename=names{i};
%         name=names{i};
%         runMStdsimple(filename,name,1,20,2^-5,fid,50);
%     end
%  end
 
% names={'AFbig','R1','R2','AFone','RX','rand','sim','AFbeta','AF225','AFpancake'};
% for i=4:7
%   filename=names{i};
%   name=names{i};
%   runMStd(filename,name,16,30,6,1);
% end

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
 