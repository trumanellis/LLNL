function hgForces=ComputeHGForces(ZoneVel)

Vx=ZoneVel(1,:);
Vy=ZoneVel(2,:);

fx=0.25*(Vx(2)-Vx(3)+Vx(4)-Vx(1));
fy=0.25*(Vy(2)-Vy(3)+Vy(4)-Vy(1));

hgForces=[
    fx,fy;
    -fx,-fy;
    fx,fy;
    -fx,-fy;]';