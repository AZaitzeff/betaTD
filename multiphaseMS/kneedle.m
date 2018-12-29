function lam=kneedle(fids,ys)
xs=fids;
N=numel(fids);
h=xs(2)-xs(1);
pp=csaps(xs,ys,1/(1 + h^3/6));
for n=1:N
    ys(n) = fnval(pp,xs(n));
end

xs=(xs-min(xs))/(max(xs)-min(xs));
ys=1-(ys-min(ys))/(max(ys)-min(ys));

yd=ys-xs;

[~,I]=max(yd);
lam=fids(I);