function options=createoptionshex(axis,odd)

if axis==-1
    options.square=1;
else
    options.square=0;
    if axis==1 && odd==0
        options.hex1f=[0,1];
        options.hex1b=[1,-1];
        options.hex2f=[1,1];
        options.hex2b=[0,-1];
        options.hex2f=[1,0];
        options.hex2b=[-1,0];
    elseif axis==1 && odd==1
        options.hex1f=[-1,1];
        options.hex1b=[0,-1];
        options.hex2f=[0,1];
        options.hex2b=[-1,-1];
        options.hex2f=[1,0];
        options.hex2b=[-1,0];
        
    elseif axis==0 && odd==0
        options.hex1f=[0,1];
        options.hex1b=[0,-1];
        options.hex2f=[1,1];
        options.hex2b=[-1,0];
        options.hex2f=[1,0];
        options.hex2b=[-1,1];
    elseif axis==0 && odd==1
        options.hex1f=[0,1];
        options.hex1b=[0,-1];
        options.hex2f=[1,0];
        options.hex2b=[-1,-1];
        options.hex2f=[1,-1];
        options.hex2b=[-1,0];
    end
end