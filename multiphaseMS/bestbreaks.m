function ind=bestbreaks(score)
best=inf;
numfids=numel(score);
ind=zeros(1,2);
for i=1:numfids-2
    for j=i+1:numfids-1
           mse=sum((score(1:i)-mean(score(1:i))).^2)+sum((score(i+1:j)-mean(score(i+1:j))).^2)+...
               sum((score(j+1:end)-mean(score(j+1:end))).^2);

       if best>mse
           best=mse;
           ind(1)=i;
           if j==(numfids-1) && (score(numfids)/score(numfids-1))>.9
               ind(2)=numfids;
           else
               ind(2)=j;
           end
       end
    end 
end