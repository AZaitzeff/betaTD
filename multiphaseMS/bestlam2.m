function [lam]=bestlam2(xs,ys)
N=numel(xs);
xtry=linspace(min(xs),max(xs),1001);
ytry=linspace(min(ys),max(ys),1001);
xstart=xs(1);
ystart=ys(1);

xend=xs(end);
yend=ys(end);
bestval=inf;
xcand=xs(2:end-1);
ycand=ys(2:end-1);
for x=xtry
    for y=ytry
        val=0;
        for n=1:(N-2)
            x0=xcand(n);
            y0=ycand(n);
            if x0<x
                x1=xstart;
                y1=ystart;
            else
                x1=xend;
                y1=yend;
            end
            f=(x1-x0)/(x1-x)*y+(x0-x)/(x1-x)*y1-y0;
            val=val+1/2*f^2;
            %valh=valh+1/2*fh^2;

        end
        if bestval>val
            bestval=val;
            lam=x;
        end

    end
end