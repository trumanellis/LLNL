function stiffM=GetStiffnessMatrix(nodes,GradBasis,quadorder)

if GradBasis(5)=='Q'
    [pts2d wgts2d]= GaussLobattoWeights2d(quadorder);

    for n=1:size(pts2d,1)
        Jmat(:,:)=JacobianMat(nodes,pts2d(n,:))';
        detJ=det(Jmat);
        invJ=inv(Jmat);

        [dwdx dwdy]=feval(GradBasis,pts2d(n,1),pts2d(n,2));
        if n>1
            stiffM=stiffM+wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        else
            stiffM=wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        end
    end
else
    [pts2d,wgts2d] = GaussLobattoWeightsTri(quadorder);
    Jmat=JacobianMatTri(nodes);
    invJ=inv(Jmat);
    detJ=det(Jmat);
    [dwdx dwdy]=GradP1Basis();

    for n=1:size(pts2d,1)
        if n>1
            stiffM=stiffM+wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        else
            stiffM=wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        end
    end
end