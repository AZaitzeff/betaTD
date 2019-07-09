function lam=highestcurvature(fids,ys)
xs=fids;
h=xs(2)-xs(1);
pp=csaps(xs,ys,1/(1 + h^3/.6));
ppd=fnder(pp,1);
ppdd=fnder(pp,2);
xtry=linspace(min(xs),max(xs),1001);
cur=zeros(1,1001);
for i=1:1001
    cur(i)=ppval(ppdd,xtry(i))/(1+ppval(ppd,xtry(i))^2);
end
plot(xtry,cur)
[~,I]=max(cur);
lam=xtry(I);
