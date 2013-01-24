function cF=ComputeGradPCornerForce(nodes,pdof,GradBasis,quadorder)

[pts2d wgts2d]= GaussLobattoWeights2d(quadorder);
if GradBasis(6)=='2'
    Dx=zeros(9,4);
    Dy=zeros(9,4);

    for n=1:size(pts2d,1)
        Jmat=JacobianMat(nodes,pts2d(n,:));
        detJ=det(Jmat);
        invJ=inv(Jmat);

        [dwdx dwdy]=feval(GradBasis,pts2d(n,1),pts2d(n,2));
        U=[dwdx dwdy];
        V=Q1dBasis(pts2d(n,1),pts2d(n,2));
        Dx=Dx+detJ*wgts2d(n)*(invJ*U')'*[1;0]*V';
        Dy=Dy+detJ*wgts2d(n)*(invJ*U')'*[0;1]*V';
    end

    cF=[Dx*pdof,Dy*pdof]';
else
    cF=zeros(2,4);

    for n=1:size(pts2d,1)
        Jmat=JacobianMat(nodes,pts2d(n,:));
        detJ=det(Jmat);
        invJ=inv(Jmat);

        [dwdx dwdy]=feval(GradBasis,pts2d(n,1),pts2d(n,2));
        cF=cF+wgts2d(n)*invJ*[dwdx dwdy]'*detJ;
    end
    cF=pdof*cF;
end