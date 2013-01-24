function Jmat=JacobianMat(nodes,localPt)

r1=nodes(1,:);
r2=nodes(2,:);
r3=nodes(3,:);
r4=nodes(4,:);
x=localPt(1);
y=localPt(2);

Jmat=[(1-y)*(r2-r1)+y*(r3-r4);
    (1-x)*(r4-r1)+x*(r3-r2)];