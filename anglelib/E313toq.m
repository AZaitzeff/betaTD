function Q=E313toq(E)
%q=q1+iq2+jq3+kq4
z=size(E,1);
Q=zeros(z,4);
Q(:,1)=cos((E(:,1)+E(:,3))/2).*cos(E(:,2)/2);
Q(:,2)=cos((E(:,1)-E(:,3))/2).*sin(E(:,2)/2);
Q(:,3)=sin((E(:,1)-E(:,3))/2).*sin(E(:,2)/2);
Q(:,4)=sin((E(:,1)+E(:,3))/2).*cos(E(:,2)/2);
