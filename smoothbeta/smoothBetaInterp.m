function u=smoothBetaInterp(guess, f, weights)

max_power = log2(size(f,3));

u = guess;

for cursize=2 .^ (3:max_power)
    f_size  = orthogonalize(resize(f, [cursize, cursize]));
%     f0_size = orthogonalize(resize(f0, [cursize, cursize]));
    u = orthogonalize(resize(u, [cursize, cursize]));
    numsteps = 0;
    
    subplot(1,2,1);
    visso3(f_size);
    subplot(1,2,2);
    ef = 0;
    ei = 1;
    
    while numsteps < 5000  && (ei - ef) > 1e-6
        [u, ei, ef] = so3implicitfid(1, 1/(cursize^2), f_size, u, 100, weights);
        numsteps = numsteps + 1;
    end
    
    % some info about each interpolation step
    fprintf('took %d steps at size %d\n', numsteps, cursize);
    fprintf('percent diff f : %f\n', sum(abs(f_size - u), 'all') / sum(abs(f_size), 'all'));
%     fprintf('percent diff f0: %f\n', sum(abs(f0_size - u), 'all') / sum(abs(f0_size), 'all'));
end

end

function q=orthogonalize(u)
    q = zeros(size(u));
    for i=1:size(u,3)
        for j=1:size(u,4)
            [Q,~] = qr(u(:,:,i,j));
            q(:,:,i,j) = Q;
        end
    end
end

function o=resize(u, newsize)
    o = zeros(3, 3, newsize(1), newsize(2));
    for i=1:3
        for j=1:3
            o(i,j,:,:) = imresize(squeeze(u(i,j,:,:)), newsize, 'bilinear');
        end
    end
end