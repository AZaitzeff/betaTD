names={'AFbig','R1','R2','AFone','RX','rand','sim','AFbeta','AF225','AFpancake','fine'};
i=1;
for fid=[200,400,600]
     filename=names{i};
     name=names{i};
     runMStdsimple(filename,name,8,16,2^-5,fid,5,200,100);
end
% % 
% i=3;
% for fid=[50]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-4,fid,30,0);
% end
% % 
% i=4;
% for fid=[150,200]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-4,fid,30,0);
% end
% 
% i=1;
% for fid=[200,250,300]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-5,fid,5,1);
% end
% % 
% % 
% i=5;
% for fid=[100]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-4,fid,30,0);
% end
% % 
% % 
% i=6;
% for fid=[25,35,50]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-4,fid,30,0);
% end
% % 
% i=7;
% for fid=[50]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,2,8,2^-4,fid,30,0);
% end

% 
% i=9;
% for fid=[100,150,200]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,4,16,2^-4,fid,40,400,0);
% end


% i=9;
% for fid=[100,150,200]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,8,16,2^-4,fid,40,1);
% end

% 
% i=10;
% for fid=[100,150]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,30,2^-4,fid,30,0);
% end


% i=2;
% for fid=[25,50,100]
%      filename=names{i};
%      name=names{i};
%      runMStdsimple(filename,name,16,20,2^-4,fid,20);
% end

%  
 
% for i=12:12
%     for fid=100:50:400
%         filename=names{i};
%         name=names{i};
%         runMStdsimple(filename,name,1,20,2^-5,fid,50);
%     end
%  end
 
% names={'AFbig','R1','R2','AFone','RX','rand','sim','AFbeta','AF225','AFpancake'};
% for i=[2:10 1]
%   filename=names{i};
%   name=names{i};
%   runMStd(filename,name,20,30,6);
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
 