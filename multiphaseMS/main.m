for enec=[.2,.1,-.1,-.2]
    for fid=[12,25,50,100,200,400]
        if enec>=0
            name=['sim' num2str(10*enec)];
        else
            name=['simn' num2str(-10*enec)];
        end
        filename='sim';
        runMStd(filename,name,fid,12,1,[8,8],1/100,1/100,2^-6,1,18,0,0,1,enec);
    end
end

for enec=[.2,.1,-.1,-.2]
    for fid=[12,25,50,100,200,400]
        if enec>=0
            name=['AFone' num2str(10*enec)];
        else
            name=['AFonen' num2str(-10*enec)];
        end
        filename='AFone';
        runMStd(filename,name,fid,12,1,[12,8],2/100,2/100,2^-6,1,18,0,0,2,enec);
    end
end

for enec=[.2,.1,0,-.1,-.2]
    for fid=[12,25,50,100,200,400]
        if enec>=0
            name=['AFbig' num2str(10*enec)];
        else
            name=['AFbign' num2str(-10*enec)];
        end
        filename='AFbig';
        runMStd(filename,name,fid,12,1,[12,8],4/100,4/100,2^-6,1,18,0,0,4,enec);
    end
end
% for dt=[.08,.16]
%    for fid=[100,200,300,400,500]
%        name=['sim' num2str(dt*100) 'dt'];
%        filename='sim';
%        runMStd(filename,name,fid,2,1,[6,6],1/100,1/100,dt,4,8,1,-1,0,1);
%    end
% end