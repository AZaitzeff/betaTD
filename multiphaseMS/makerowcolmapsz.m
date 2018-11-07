function [xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,rl,cl,m,n)
xdir = zeros(4000,801);
xsizes=zeros(4000,1);
step=1;
for i=1:m
    [arrays,count,sizes]=decomposez(mask(i,:));
    for z=1:count
        xdir(step,1)=i;
        p=sizes(z);
        xdir(step,2:p+1)=arrays(z,1:p);
        xsizes(step)=p;
        step=step+1;
        if step>4000
            'warning too large'
            break
            
        end
    end
    if step>4000
        break
    end
end
ydir = zeros(4000,801);
ysizes=zeros(4000,1);
step=1;
for j=1:n
    [arrays,count,sizes]=decomposez(mask(:,j)');
    for z=1:count
        ydir(step,1)=j;
        p=sizes(z);
        ydir(step,2:p+1)=arrays(z,1:p);
        ysizes(step)=p;
        step=step+1;
        if step>4000
            break
        end
    end
    if step>4000
        'warning too large'
        break
    end
end