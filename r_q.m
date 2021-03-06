function q=r_q(R)

%   q=r_q(R) translates the 3x3 DCM R into the corresponding quaternion q

tr=trace(R);
b2(1) = sqrt((1+tr)/4);
b2(2) = sqrt((1+2*R(1,1)-tr)/4);
b2(3) = sqrt((1+2*R(2,2)-tr)/4);
b2(4) = sqrt((1+2*R(3,3)-tr)/4);

[v,i] = max(b2);
switch i,
case 1
    b(1) = (b2(1));
    b(2) = (R(2,3) - R(3,2))/4/b(1);
    b(3) = (R(3,1) - R(1,3))/4/b(1);
    b(4) = (R(1,2) - R(2,1))/4/b(1);
case 2
    b(2) = (b2(2));
    b(1) = (R(2,3) - R(3,2))/4/b(2);
    if b(1)<0, b(2)=-b(2); b(1)=-b(1); end
    b(3) = (R(1,2) + R(2,1))/4/b(2);
    b(4) = (R(3,1) + R(1,3))/4/b(2);
case 3
    b(3) = (b2(3));
    b(1) = (R(3,1) - R(1,3))/4/b(3);
    if b(1)<0, b(3)=-b(3); b(1)=-b(1); end
    b(2) = (R(1,2) + R(2,1))/4/b(3);
    b(4) = (R(2,3) + R(3,2))/4/b(3);
case 4
    b(4) = (b2(4));
    b(1) = (R(1,2) - R(2,1))/(4*b(4));
    if b(1)<0, b(4)=-b(4); b(1)=-b(1); end
    b(2) = (R(3,1) + R(1,3))/(4*b(4));
    b(3) = (R(2,3) + R(3,2))/(4*b(4));
end

q = [b(2:4) b(1)];

