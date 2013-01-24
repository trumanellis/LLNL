%% Compute the hydro corner forces
% These are simply the values of the nodal forces at each node of a
% given zone
if ~isPresRefined
    for i=1:NZ
        cornerForce(:,:,i)=ComputeGradPCornerForce(oldmesh(:,:,i),pressureZ(i),GradBasis,stiffQuadOrder);
    end
else
    for i=1:NZ
        cornerForce(:,:,i)=ComputeGradPCornerForce(oldmesh(:,:,i),pressureQ1d(:,i),GradBasis,stiffQuadOrder);
    end
end

if hgfrac > 0 && strcmp(Method,'Q1Q0')
    % Compute the nodal hour glass forces
    % These are the artificial forces that are designed to "resist" the
    % hourglass modes
    % The magnitude of this force can be scaled with the parameter hgfrac
    for i=1:NZ
        zonevel=OLDvelocityQuad(:,topo(:,i));

        thisHGForce=1/(4*2*dt)*massZ(i)*ComputeHGForces(zonevel);

        % Set the hourglass corner force
        cornerForce(:,:,i)=cornerForce(:,:,i)+hgfrac*thisHGForce(:,:);
    end
end

% These are the tensor artificial viscosity forces
if Qfrac ==1
    for i=1:NZ
        zonevel=OLDvelocityQuad(:,Quadmap(:,i));
        if ~isPresRefined
            % Calculate the divergence
            dV=LocalInterpolateVelocityGrad(zonevel,oldmesh(:,:,i),[0.5,0.5],GradBasis);
            DivV_avg=dV(1,1)+dV(2,2);
            if DivV_avg<0
                % The viscosity coefficient
                mu=OLDdensityZ(i)*len(i)*(-qquad*DivV_avg*len(i)+qlin*soundSpeedZ(i));
                % Calculate a local stiffness matrix
                stiffmat=GetStiffnessMatrix(oldmesh(:,:,i),GradBasis,stiffQuadOrder);
                thisTenQForce=mu*(stiffmat*zonevel');
                cornerForce(:,:,i)=cornerForce(:,:,i)-thisTenQForce';
            end
        else
            thisTenQForce=ComputeTenQForce(oldmesh(:,:,i),zonevel,OLDdensityQ1d(:,i),pressureQ1d(:,i),gamma,len(i),GradBasis,qquad,qlin,stiffQuadOrder);
            cornerForce(:,:,i)=cornerForce(:,:,i)-thisTenQForce;
        end
    end
end

% Sum-in hydro corner forces to obtain the hydro nodal forces
forceN=zeros(2,ndofpD);
for i=1:NZ
    forceN(:,Quadmap(:,i))=forceN(:,Quadmap(:,i))+cornerForce(:,:,i);
end