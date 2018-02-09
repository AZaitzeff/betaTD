function Q=E313toq(E)
%q=iq1+jq2+kq3+q4
z=size(E,1);
Q=zeros(z,4);
Q(:,1)=cos((E(:,1)-E(:,3))/2).*sin(E(:,2)/2);
Q(:,2)=sin((E(:,1)-E(:,3))/2).*sin(E(:,2)/2);
Q(:,3)=sin((E(:,1)+E(:,3))/2).*cos(E(:,2)/2);
Q(:,4)=cos((E(:,1)+E(:,3))/2).*cos(E(:,2)/2);