function [xdir,ydir]=makerowcolmaps(mask,rl,cl,m,n)

xdir = cell(1,rl);
step=1;
for i=1:m
    [arrays,count]=decompose(mask(i,:));
    if count>0
        xdir{step}=[i,count,arrays];
        step=step+1;
    end
end
ydir = cell(1,cl);
step=1;
for j=1:n
    [arrays,count]=decompose(mask(:,j)');
    if count>0
        ydir{step}=[j,count,arrays];
        step=step+1;
    end
end