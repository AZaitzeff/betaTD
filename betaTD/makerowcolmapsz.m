function [xdir,ydir,xsizes,ysizes]=makerowcolmapsz(mask,m,n)
sizex=[8*m,n+1];
xdir = zeros(sizex);
xsizes=zeros(sizex(1),1);
coder.varsize('xdir',[],[1,1]);
coder.varsize('xsizes',[],[1 0]);
step=1;
for i=1:m
    [arrays,count,sizes]=decomposez(mask(i,:),sizex(2));
    maxp=0;
    for q=1:count
        maxp=max(sizes(q),maxp);
    end
    fstep=step+count;
    if maxp>=sizex(2) && sizex(1)<fstep
        sizex(1)=sizex(1)*2;
        sizex(2)=ceil(maxp*1.1);
        [xdir,xsizes]=growarray(xdir,xsizes,sizex(1)/2,sizex,0,1);
    elseif maxp>=sizex(2)
        sizex(2)=ceil(maxp*1.1);
        [xdir,xsizes]=growarray(xdir,xsizes,sizex(1),sizex,0,1);
    elseif sizex(1)<fstep
        sizex(1)=sizex(1)*2;
        [xdir,xsizes]=growarray(xdir,xsizes,sizex(1)/2,sizex,0,1);
    end
    for z=1:count
        xdir(step,1)=i;
        p=sizes(z);
        
        
        xdir(step,2:p+1)=arrays(z,1:p);
        xsizes(step)=p;
        step=step+1;
    end
end
sizey=[8*n,m+1];
ydir = zeros(sizey);
ysizes=zeros(sizey(1),1);
coder.varsize('ydir',[],[1,1]);
coder.varsize('ysizes',[],[1 0]);
step=1;
for j=1:n
    [arrays,count,sizes]=decomposez(mask(:,j),sizey(2));
    maxp=0;
    for q=1:count
        maxp=max(sizes(q),maxp);
    end
    fstep=step+count;
    if maxp>=sizey(2) && sizey(1)<fstep
        sizey(1)=sizey(1)*2;
        sizey(2)=ceil(maxp*1.1);
        [ydir,ysizes]=growarray(ydir,ysizes,sizey(1)/2,sizey,0,1);
    elseif maxp>=sizey(2)
        sizey(2)=ceil(maxp*1.1);
        [ydir,ysizes]=growarray(ydir,ysizes,sizey(1),sizey,0,1);
    elseif sizey(1)<fstep
        sizey(1)=sizey(1)*2;
        [ydir,ysizes]=growarray(ydir,ysizes,sizey(1)/2,sizey,0,1);
    end
    
    for z=1:count
        ydir(step,1)=j;
        p=sizes(z);
        ydir(step,2:p+1)=arrays(z,1:p);
        ysizes(step)=p;
        step=step+1;
    end
end