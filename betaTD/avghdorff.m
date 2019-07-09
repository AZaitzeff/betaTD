function val=avghdorff(first,second,M,N)
n1=size(first,1);
n2=size(second,1);
val1=0;
for p=first
    val1=val1+min([sqrt(min((p-second).^2)),p(1),M-p(1),p(2),N-p(2)]);
end

val2=0;

for p=second
    val2=val2+min([sqrt(min((p-first).^2)),p(1),M-p(1),p(2),N-p(2)]);
end

val=val1/n1+val2/n2;