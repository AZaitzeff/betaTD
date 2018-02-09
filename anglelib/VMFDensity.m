function [y] = VMFDensity(x, mu, k)
    mu = mu(:)';
    p=size(mu,2);
    if(size(x,2)~=p)
        error('the dimension of sample and mean should be consistent');
    end
%     y = Cp(k,p)*exp(k*(mu*x')');
    
    y = exp(logCp(k,p)+k*(mu*x')');
end