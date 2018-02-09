function Q=randgridquat(N)
    Q=randn(4,N);
    normQ=vecnorm(Q);
    Q=Q'./(normQ'*[1,1,1,1]);
end