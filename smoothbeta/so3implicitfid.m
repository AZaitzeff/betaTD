function [u,ei,ef]=so3implicitfid(nt,dt,f,uin,lam)
    ei = 0;
    ef = 0;

    m = size(f,3);
    n = size(f,4);

    I=sqrt(-1);
    wx=exp(I*2*pi/n);       % nth root of unity.
    wy=exp(I*2*pi/m);       % mth root of unity.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient matrix for taking -laplacian:
    % M = f;
    [x,y]=meshgrid(1:n,1:m);
%     x=x';
%     y=y';
    M = 4-wx.^(x-1)-wx.^(1-x)-wy.^(y-1)-wy.^(1-y);

    M=M*max(n^2,m^2);        % Divide by dx^2.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     u = zeros(3,3,m,n);
%     for i = 1:m
%         for j = 1:n
%             u(:,:,i,j) = eye(3);
%         end
%     end
    u = uin;
    
    for t = 1:nt
        % calculate lagrange multiplier
        [lagrange,e] = calclagmult(u, f, lam);
        rhs = u - dt * pixelmult(u, lagrange);
        
        if t == 1
            ei = e;
        end
        
        % fidelity term
        rhs = rhs + dt * lam * (f - u);
        
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
    
    [~, ef] = calclagmult(u, f, lam);
end
