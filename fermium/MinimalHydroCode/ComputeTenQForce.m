function thisTenQForce=ComputeTenQForce(nodes,zonevel,subdensity,subpressure,gamma,len,GradBasis,qquad,qlin,quadorder)

if GradBasis(6)=='1'
    len=1/sqrt(ComputeInvLengthSquared(nodes));
    dV=LocalInterpolateVelocityGrad(zonevel,nodes,[.5 .5],GradBasis);
    DivV=dV(1,1)+dV(2,2);
    coef=subdensity.*len.*(-qquad*DivV.*len+qlin*sqrt(gamma*subpressure./subdensity));
    
    [pts2d wgts2d]= GaussLobattoWeights2d(quadorder);
    for n=1:size(pts2d,1)
        Jmat(:,:)=JacobianMat(nodes,pts2d(n,:))';
        detJ=det(Jmat);
        invJ=inv(Jmat);

        [dwdx dwdy]=feval(GradBasis,pts2d(n,1),pts2d(n,2));
        if n>1
            stiffM=stiffM+coef*wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        else
            stiffM=coef*wgts2d(n)*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
        end
    end
    thisTenQForce=zonevel*stiffM;
else
%     quadnodes=zeros(2,4,4);
%     quadnodes(:,:,1)=[nodes(:,1),...
%         (nodes(:,2)+nodes(:,1))/2,...
%         mean(nodes,2),...
%         (nodes(:,4)+nodes(:,1))/2];
%     quadnodes(:,:,2)=[(nodes(:,2)+nodes(:,1))/2,...
%         nodes(:,2),...
%         (nodes(:,3)+nodes(:,2))/2,...
%         mean(nodes,2)];
%     quadnodes(:,:,3)=[mean(nodes,2),...
%         (nodes(:,3)+nodes(:,2))/2,...
%         nodes(:,3),...
%         (nodes(:,3)+nodes(:,4))/2,...
%         ];
%     quadnodes(:,:,4)=[(nodes(:,4)+nodes(:,1))/2,...
%         mean(nodes,2),...
%         (nodes(:,3)+nodes(:,4))/2,...
%         nodes(:,4)];

%     len=zeros(4,1);
%     for i=1:4
%         len(i)=1/sqrt(ComputeInvLengthSquared(quadnodes(:,:,i)));
%     end
    
%     dV=zeros(2,2,4);
%     dV(:,:,1)=LocalInterpolateVelocityGrad(zonevel,nodes,[.25 .25],GradBasis);
%     dV(:,:,2)=LocalInterpolateVelocityGrad(zonevel,nodes,[.75 .25],GradBasis);
%     dV(:,:,3)=LocalInterpolateVelocityGrad(zonevel,nodes,[.75 .75],GradBasis);
%     dV(:,:,4)=LocalInterpolateVelocityGrad(zonevel,nodes,[.25 .75],GradBasis);
    
%     DivV=dV(1,1,:)+dV(2,2,:);
%     DivV=DivV(:);

%     mu=subdensity.*len.*(-qquad*DivV.*len+qlin*sqrt(gamma*subpressure./subdensity));
    len=0.5*len;
    [pts2d wgts2d]= GaussLobattoWeights2d(quadorder);

    stiffM=zeros(9,9);
    
    for n=1:size(pts2d,1)
        dV=LocalInterpolateVelocityGrad(zonevel,nodes,pts2d(n,:),GradBasis);
        DivV=dV(1,1)+dV(2,2);
        wQ1d=Q1dBasis(pts2d(n,1),pts2d(n,2));
        mu=(subdensity'*wQ1d)*len*(-qquad*DivV*len+qlin*sqrt(gamma*(subpressure'*wQ1d)/(subdensity'*wQ1d)));
        
        Jmat=JacobianMat(nodes,pts2d(n,:))';
        detJ=det(Jmat);
        invJ=inv(Jmat);

        [dwdx dwdy]=feval(GradBasis,pts2d(n,1),pts2d(n,2));
        stiffM=stiffM+wgts2d(n)*mu*([dwdx dwdy]*invJ)*([dwdx dwdy]*invJ)'*detJ;
    end
    thisTenQForce=zonevel*stiffM;
end