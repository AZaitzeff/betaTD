function Pm=getsymmetries(s1)
    %cubic
    if strcmp( s1,'cube') || strcmp( s1,'cubic')
        sym = CubSymmetries();
    %Hexigon
    else
        sym = HexSymmetries();
    end
    N=size(sym,2);
    Pm=zeros(4,4,N);
    for i=1:N
        Pm(:,:,i)=quatmatrix2(sym(:,i));
    end
end