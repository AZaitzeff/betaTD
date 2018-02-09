function E=qtoE313(Q)
%q=iq1+jq2+kq3+q4
z=size(Q,1);
E=zeros(z,3);
E(:,1)=atan2(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4),-(Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)));
E(:,2)=acos(-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2+Q(:,4).^2);
E(:,3)=atan2(Q(:,1).*Q(:,3)-Q(:,2).*Q(:,4),Q(:,2).*Q(:,3)+Q(:,1).*Q(:,4));

E(sign(E)<0) = E(sign(E)<0) + 2*pi;