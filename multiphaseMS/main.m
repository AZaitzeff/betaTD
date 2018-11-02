for dt=[.08,.16]
    for fid=[250,300,350]
        name=['AFone' num2str(dt*100) 'dt' num2str(fid) 'fid'];
        filename='AFone';
        runMStd(filename,name,fid,12,1,[11,7],1/100,1/100,dt,4,12,1,1);
    end
end

for step=[4,8,16]
for dt=[.08,.16,.32]
    for fid=[250,300,350]
        name=['AFbig' num2str(dt*100) 'dt' num2str(fid) 'fid' num2str(step) 'step'];
        filename='AFbig';
        runMStd(filename,name,fid,12,1,[12,7],1/100,1/100,dt,step,12,1,0);
    end
end
end


%filename='sim';
%runMStd(filename,'simseg',[50,100,150,200],1000,6,1,Ks,.02);