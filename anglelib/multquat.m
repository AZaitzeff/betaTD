function  n = multquat(r,q)
    Q=[[r(4),r(3),-r(2),r(1)];[-r(3),r(4),r(1),r(2)];[r(2),-r(1),r(4),r(3)];[-r(1),-r(2),-r(3),r(4)]];
    n=Q*q;
end