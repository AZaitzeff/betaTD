function Um=getalpha2alpha()
    %cubic
    uniquemis = [
     1.0000    0.9958    0.7071    0.8514    0.8624    0.8660
          0    0.0000    0.2959   -0.5000   -0.4979    0.2500
          0    0.0000   -0.6422    0.0000    0.0459    0.4330
          0   -0.0918   -0.0000   -0.1582   -0.0795   -0.0000];
    N=size(uniquemis,2);
    Um=zeros(4,4,N);
    for i=1:N
        Um(:,:,i)=quatmatrix2(uniquemis(:,i));
    end
end
