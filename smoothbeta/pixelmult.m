function c=pixelmult(a, b)
    c = zeros(size(a));
    assert(isequal(size(a), size(b)));

    [~, ~, m, n] = size(a);
    for i = 1:m
        for j = 1:n
            c(:,:,i,j) = a(:,:,i,j) * b(:,:,i,j);
        end
    end
end
