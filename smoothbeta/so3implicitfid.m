function [u,ei,ef]=so3implicitfid(nt,dt,f,uin,lam, w)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computes the smooth approximation of f using gradient descent
    % Arguments:
    %   nt  - integer - number of time steps
    %   dt  - double  - size of each time step
    %   f   - 3 x 3 x n x m double array - SO(3) valued target "image"
    %   uin - 3 x 3 x n x m double array - SO(3) valued initial guess
    %   lam - double  - non spatially varying fidelity weight
    %   w   - n x m double array - spatially varying fidelity weight
    %           (optional)
    %
    % Note: lam and w eventually just get multiplied together, this is just
    %     a holdover from before when there was no spatially varying weight
    %
    % Returns: u - 3 x 3 x n x m double array - SO(3) valued result of
    %     gradient descent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin > 5
        weights = w;
    else
        weights = ones(size(f, 3), size(f, 4));
    end
    
    ei = 0;

    n = size(f,3);
    m = size(f,4);

    I=sqrt(-1);
    wx=exp(I*2*pi/n);       % nth root of unity.
    wy=exp(I*2*pi/m);       % mth root of unity.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient matrix for taking -laplacian:
    % M = f;
    [x,y]=meshgrid(1:n,1:m);
    x=x';
    y=y';
    M = 4-wx.^(x-1)-wx.^(1-x)-wy.^(y-1)-wy.^(1-y);

    M=M*max(n^2,m^2);        % Divide by dx^2.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = uin;
    
    for t = 1:nt
        % calculate lagrange multiplier
        [lagrange,e] = calclagmult(u, f, lam, weights);
        rhs = u - dt * pixelmult(u, lagrange);
        
        if t == 1
            ei = e;
        end
        
        % fidelity term
        fid = (f - u);
        for i=1:n
            for j=1:m
                fid(:,:,i,j) = fid(:,:,i,j) * weights(i,j);
            end
        end
        rhs = rhs + dt * lam * fid;
        
        % laplacian term
        for i=1:3
            for j=1:3
                rhsf = fft2(squeeze(rhs(i,j,:,:)));
                rhsf = rhsf ./ (1 + dt * M);
                rhs(i,j,:,:) = real(ifft2(rhsf));
            end
        end
        
        % debugging print statments and showing the figure
%         if mod(t,10) == 0
%             visso3(rhs);
%             saveas(gcf, strcat('../pics/', num2str(t, '%05d'), '.png'));
%         end
%         fprintf('step: %d \tenergy: %e \torthogs: %d \tsum: %e\n',...
%             t, e, countorthogs(rhs), sum(rhs .^ 2, 'all'));
        if isnan(e)
            error("blew up");
        end
        
        % update and repeat
        u = rhs;
    end
    
    [~, ef] = calclagmult(u, f, lam, weights);
end
