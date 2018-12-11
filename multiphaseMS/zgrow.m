function [newrow,newcol,stop]=zgrow(row,col,width,M,N)
xneigh=[1,0,-1,0];
yneigh=[0,1,0,-1];
newrow=zeros(1,M*N);
newcol=zeros(1,M*N);
z=numel(row);
newrow(1:z)=row;
newcol(1:z)=col;
W=ones(M+2,N+2,'logical');
W(:,1)=0;
W(1,:)=0;
W(:,N+2)=0;
W(M+2,:)=0;
for i=1:z
    W(row(i)+1,col(i)+1)=0;
end
start=1;
stop=z;
z=z+1;
for D=1:width
    for i=start:stop
        y=newrow(i);
        x=newcol(i);
        for j=1:4
            xc=x+xneigh(j);
            yc=y+yneigh(j);
            if W(yc+1,xc+1)
                W(yc+1,xc+1)=0;
                newrow(z)=yc;
                newcol(z)=xc;
                z=z+1;
            end
            
        end    
    end

    start=stop+1;
    stop=z-1;
end

newrow=newrow(1:stop);
newcol=newcol(1:stop);
end