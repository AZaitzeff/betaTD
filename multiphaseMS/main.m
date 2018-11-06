%  for dt=[.08,.16]
%      for fid=[250,300,350]
%          name=['AFone' num2str(dt*100) 'dt'];
%          filename='AFone';
%          runMStd(filename,name,fid,12,1,[11,7],1/100,1/100,dt,4,24,1,-1,0,1);
%      end
% end
for step=[4]
for dt=[.16]
    for fid=[300,400]
        name=['AFbig' num2str(dt*100) 'dt' num2str(step) 'step'];
        filename='AFbig';
        runMStd(filename,name,fid,12,1,[6,10],1/100,1/100,dt,step,18,1,.5,2,1);
    end
end
end
for step=[4]
for dt=[.16]
    for fid=[300,400]
        name=['AFbig' num2str(dt*100) 'dt' num2str(step) 'stepl'];
        filename='AFbig';
        runMStd(filename,name,fid,12,1,[6,10],1/100,1/100,dt,step,18,1,.5,-1,1);
    end
end
end
% for dt=[.08,.16]
%    for fid=[100,200,300,400,500]
%        name=['sim' num2str(dt*100) 'dt'];
%        filename='sim';
%        runMStd(filename,name,fid,2,1,[6,6],1/100,1/100,dt,4,8,1,-1,0,1);
%    end
% end