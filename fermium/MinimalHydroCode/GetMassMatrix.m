function massM=GetMassMatrix(nodes,Basis,order,densityQ1d)

[pts2d,wgts2d] = GaussLobattoWeights2d(order);

if nargin < 4
    for n=1:size(pts2d,1)
        Jmat=JacobianMat(nodes,pts2d(n,:));
        detJ=det(Jmat);
        w=feval(Basis,pts2d(n,1),pts2d(n,2));
        if n>1
            massM=massM+wgts2d(n)*(w*w')*detJ;
        else
            massM=wgts2d(n)*(w*w')*detJ;
        end
    end
else
    massM=zeros(9);
    for n=1:size(pts2d,1)
        Jmat=JacobianMat(nodes,pts2d(n,:));
        detJ=det(Jmat);
        w=feval(Basis,pts2d(n,1),pts2d(n,2));
        wQ1d=Q1dBasis(pts2d(n,1),pts2d(n,2));
        
        massM=massM+wgts2d(n)*(densityQ1d'*wQ1d)*(w*w')*detJ;
    end
end