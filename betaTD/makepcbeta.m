n=size(dict,2);

pcBeta=zeros(400,400,3);
temp1=zeros(400,400);
temp2=zeros(400,400);
temp3=zeros(400,400);
for z=1:n
    temp=qtoE313(dict(:,z));
    temp1(mapall==z)=temp(1);
    temp2(mapall==z)=temp(2);
    temp3(mapall==z)=temp(3);
end

pcBeta(:,:,1)=temp1;
pcBeta(:,:,2)=temp2;
pcBeta(:,:,3)=temp3;