%  for dt=[.08,.16]
%      for fid=[250,300,350]
%          name=['AFone' num2str(dt*100) 'dt'];
%          filename='AFone';
%          runMStd(filename,name,fid,12,1,[11,7],1/100,1/100,dt,4,24,1,-1,0,1);
%      end
% end
for Ks=[[4,8];[5,9];[6,10]]'
for step=[4]
for dt=[.08,.16]
    for fid=[150,200,250]
        name=['AFbig' num2str(dt*100) 'dt' num2str(step) 'step' num2str(Ks(1)) 'Ks'];
        filename='AFbig';
        runMStd(filename,name,fid,12,1,Ks,1/100,1/100,dt,step,24,1,.5,2,1);
    end
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