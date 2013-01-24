function dV=LocalInterpolateVelocityGrad(dof,nodes,localPoint,GradBasis)

Jmat=JacobianMat(nodes,localPoint);
Jinv=inv(Jmat);
[dwdx dwdy]=feval(GradBasis,localPoint(1),localPoint(2));

interpMat=Jinv*[dwdx dwdy]';
dV=dof*interpMat';