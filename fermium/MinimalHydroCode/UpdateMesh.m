% UpdateMesh
% Integrate accelerations to get velocities
NEWvelocityQuad=OLDvelocityQuad+dt*accelerationQuad;

newnodes=oldnodes+dt*NEWvelocityQuad(:,1:NN);
newmesh=ComputeMesh(newnodes,topo);

if strcmp(Method,'Q2Q1d')
    if isEdgeCenIndependent == true
        % This block moves the edge and center nodes
        newedgenodes=oldedgenodes+dt*NEWvelocityQuad(:,NN+1:NN+NE);
        newcennodes=oldcennodes+dt*NEWvelocityQuad(:,NN+NE+1:end);

        % Project new edge locations back onto line between corner nodes
        % (Preserves Quad shape)
        for i=1:NE
            p1=newnodes(:,edgelist(2,i));
            p0=newnodes(:,edgelist(1,i));
            q=newedgenodes(:,i);
            newedgenodes(:,i)=[p1(1)-p0(1), p1(2)-p0(2); p0(2)-p1(2), p1(1)-p0(1)]\[q(1)*(p1(1)-p0(1))+q(2)*(p1(2)-p0(2)); p0(2)*(p1(1)-p0(1))-p0(1)*(p1(2)-p0(2))];
        end
    else
        % This block just averages the locations of the edge and center nodes
        newedgenodes=(newnodes(:,edgelist(1,:))+newnodes(:,edgelist(2,:)))/2;
        for i=1:NZ
            newcennodes(:,i)=ComputeCentroid(newmesh(:,:,i));
        end
    end

    allnodes=[newnodes,newedgenodes,newcennodes];
end