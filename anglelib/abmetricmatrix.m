function dist=abmetricmatrix(alpha,beta)
    transfmatrices();
    

    dist=inf;
    for i=1:6
        for j=1:24
            na=T{i}*alpha;
            nb=Pmb{j}*beta;
            g=nb/na;
            dist=min(dist,acos((trace(g)-1)/2));
        end
    end
end