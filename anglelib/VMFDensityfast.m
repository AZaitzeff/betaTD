function [y] = VMFDensityfast(x, mu, k)
    mu = mu(:)';
    p=size(mu,2);
    if(size(x,2)~=p)
        error('the dimension of sample and mean should be consistent');
    end
%     y = Cp(k,p)*exp(k*(mu*x')');
    if k<700
        val = log((2*pi)^(-p/2) * (k.^(p/2-1) ./ besseli(p/2-1, k)));
    else
        y2=log((2*pi)^(-p/2) * (700.^(p/2-1) ./ besseli(p/2-1, 700)));
        y1=log((2*pi)^(-p/2) * (699.^(p/2-1) ./ besseli(p/2-1, 699)));
        val = y2+(y2-y1)*(k-700);
    end
    y = real(exp(val+k*(mu*x')'));
end