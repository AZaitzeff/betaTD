function ind=bestbreaks(score)
best=inf;
copy=score;
numfids=numel(score);

score(1)=(copy(1)+copy(2))/2;
score(end)=(copy(end)+copy(end))/2;
for n=2:numfids-1
   score(n)=median(copy(n-1:n+1)); 
end
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