function minvec=bruteforceavg(alpha,Q)
    n=size(alpha,1);
    %score=alpbmetric(alpha,Q);
    score=abmetric2(alpha,Q);
    score=sum(score,2)/n;
    [~,I]=min(score);
    minvec=Q(I,:);
end
    
    