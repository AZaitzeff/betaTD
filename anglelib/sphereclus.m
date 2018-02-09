function [mu]=sphereclus(x,k)
    xAp=0.001:0.1:700;
    yAp = Ap(xAp);
    
    mu0=randgridquat(k);
    kappa=50*ones(1,k);
    alpha=ones(1,k)/k;
    T=alphatobetatrans();
    Pm=getsymmetries('cubic');
    [n,~]=size(x);
    m=144;
    r=zeros(n,m,k);
    maxiter=10000;
    tol=1e-6;
    for nt=1:maxiter
        mu=mu0;
        for c=1:k
            for z=1:m
                [i,j]=num2ind(z);
                r(:,z,k)=alpha(c)*phi(x',T(:,:,i)'*Pm(:,:,j)*mu(k,:)',kappa(c));
            end
        end
        r = bsxfun(@rdivide,r,sum(sum(r,3),2));
        gamma=zeros(k,4);
        for c=1:k
            for z=1:m
                [i,j]=num2ind(z);
                gamma(c,:)=gamma(c,:)+sum(bsxfun(@times,(Pm(:,:,j)'*T(:,:,i)*x')',r(:,z,c)),1);
            end
            alpha(c)=sum(sum(r(:,:,c)));
            gam2=norm(gamma(c,:));
            mu(c,:)=gamma/gam2;
            kappa(c) = interp1(yAp, xAp, (gam2/(alpha(c)*n)), 'linear', 'extrap');

        end
        sum(abs(mu0(:)-mu(:)))
        if sum(abs(mu0(:)-mu(:)))<tol
            break
        end
    end
    %clusters on the 3 sphere k groups
    
end

function val=phi(x,mu,kappa)
    const=sqrt(kappa)/((2*pi)^(3/2)*besseli(1/2,kappa));
    val=const*exp(kappa*mu'*x);
end

function val=Ap(u)
    val=besseli(3/2,u)./besseli(1/2,u);
end

function [i,j]=num2ind(z)
    i=mod(z-1,6)+1;
    j=ceil(z/6);
end