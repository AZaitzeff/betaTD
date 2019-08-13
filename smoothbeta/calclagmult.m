function [g,e]=calclagmult(u, f, lam, weights)
    % Computes the lagrange multiplier, g.
    % (Also computes the current value of the energy function using the forward
    % difference, but this is just out of convenience)

    ub = buffer(u);
    [fdx, bdx] = xdiff(ub);
    [fdy, bdy] = ydiff(ub);

    e = sum(fdx(:) .^ 2 + lam * (f(:) - u(:)) .^ 2); 
    g = -0.5 * (transpixelmult(fdx) + transpixelmult(bdx) +...
        transpixelmult(fdy) + transpixelmult(bdy));
    
    g = g + fidelcorrection(u, f, lam, weights);
end

function o=fidelcorrection(u, f, lam, weights)
    assert(isequal(size(u),size(f)));
    o = zeros(size(u));
    if lam == 0
        return
    end
    
    n = size(u,3);
    m = size(u,4);
    for i = 1:n
        for j = 1:m
            o(:,:,i,j) = (u(:,:,i,j)' * f(:,:,i,j) ... 
                          + f(:,:,i,j)' * u(:,:,i,j)) * 0.5 ...
                         - eye(3);
            o(:,:,i,j) = o(:,:,i,j) * weights(i,j);
        end
    end
    o = o * lam;
end

function c=transpixelmult(a)
    c = zeros(size(a));

    n = size(a,3);
    m = size(a,4);
    for i = 1:n
        for j = 1:m
            c(:,:,i,j) = a(:,:,i,j)' * a(:,:,i,j);
        end
    end
end

function ub=buffer(u)
    [c,d,n,m]=size(u);

    ub=zeros(c,d,n+2,m+2);	% Initialize ub.

    % Fill in the buffer:
    ub(:,:,2:n+1,2:m+1)=u;	% Bulk part of ub is just u.
    ub(:,:,2:n+1,1)=u(:,:,1:n,m);
    ub(:,:,2:n+1,m+2)=u(:,:,1:n,1);
    ub(:,:,1,2:m+1)=u(:,:,n,1:m);
    ub(:,:,n+2,2:m+1)=u(:,:,1,1:m);

    ub(:,:,1,1) = u(:,:,n,m);
    ub(:,:,1,m+2) = u(:,:,n,1);
    ub(:,:,n+2,1) = u(:,:,1,m);
    ub(:,:,n+2,m+2) = u(:,:,1,1);
end

function [fdx, bdx]=xdiff(ub)
    % Assumes ub is a buffered matrix-valued image with a buffer of size 1
    % Calculates either forward or backward difference quotient in the
    % x direction depending on the value of forward
    n = size(ub,3) - 2; % Subtract off buffer size
    m = size(ub,4) - 2;
    fdx = (ub(:,:,2:n+1,3:m+2) - ub(:,:,2:n+1,2:m+1)) * m;
    bdx = (ub(:,:,2:n+1,2:m+1) - ub(:,:,2:n+1,1:m)) * m;
end

function [fdy, bdy]=ydiff(ub)
    % Assumes ub is a buffered matrix-valued image with a buffer of size 1
    % Calculates either forward or backward difference quotient in the
    % y direction depending on the value of forward
    n = size(ub,3) - 2; % Subtract off buffer size
    m = size(ub,4) - 2;
    fdy = (ub(:,:,3:n+2,2:m+1) - ub(:,:,2:n+1,2:m+1)) * n;
    bdy = (ub(:,:,2:n+1,2:m+1) - ub(:,:,1:n,2:m+1)) * n;
end
