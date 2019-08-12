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
    
    m = size(u,3);
    n = size(u,4);
    for i = 1:m
        for j = 1:n
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

    m = size(a,3);
    n = size(a,4);
    for i = 1:m
        for j = 1:n
            c(:,:,i,j) = a(:,:,i,j)' * a(:,:,i,j);
        end
    end
end

function ub=buffer(u)
    [c,d,m,n]=size(u);

    ub=zeros(c,d,m+2,n+2);	% Initialize ub.

    % Fill in the buffer:
    ub(:,:,2:m+1,2:n+1)=u;	% Bulk part of ub is just u.
    ub(:,:,2:m+1,1)=u(:,:,1:m,n);
    ub(:,:,2:m+1,n+2)=u(:,:,1:m,1);
    ub(:,:,1,2:n+1)=u(:,:,m,1:n);
    ub(:,:,m+2,2:n+1)=u(:,:,1,1:n);

    ub(:,:,1,1) = u(:,:,m,n);
    ub(:,:,1,n+2) = u(:,:,m,1);
    ub(:,:,m+2,1) = u(:,:,1,n);
    ub(:,:,m+2,n+2) = u(:,:,1,1);
end

function [fdx, bdx]=xdiff(ub)
    % Assumes ub is a buffered matrix-valued image with a buffer of size 1
    % Calculates either forward or backward difference quotient in the
    % x direction depending on the value of forward
    m = size(ub,3) - 2; % Subtract off buffer size
    n = size(ub,4) - 2;
    fdx = (ub(:,:,2:m+1,3:n+2) - ub(:,:,2:m+1,2:n+1)) * n;
    bdx = (ub(:,:,2:m+1,2:n+1) - ub(:,:,2:m+1,1:n)) * n;
end

function [fdy, bdy]=ydiff(ub)
    % Assumes ub is a buffered matrix-valued image with a buffer of size 1
    % Calculates either forward or backward difference quotient in the
    % y direction depending on the value of forward
    m = size(ub,3) - 2; % Subtract off buffer size
    n = size(ub,4) - 2;
    fdy = (ub(:,:,3:m+2,2:n+1) - ub(:,:,2:m+1,2:n+1)) * m;
    bdy = (ub(:,:,2:m+1,2:n+1) - ub(:,:,1:m,2:n+1)) * m;
end
