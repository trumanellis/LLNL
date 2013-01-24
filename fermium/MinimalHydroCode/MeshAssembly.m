%% Mesh Assembly
%Compute the mesh nodes
refnodes=ComputeReferenceMeshNodes(NZx,NZy,[xmin,xmax],[ymin,ymax]);

%Compute the mesh connectivity
topo=ComputeMeshTopology(NZx,NZy);

% =============================== PERTURBED MESH ===========================

if jitter > 0
    %Construct the node weights: i.e. the numer of zones which share each node
    nodeweights=ComputeNodeWeights(topo);

    %Compute the boundary and interior nodes
    [bnodes, interiornodes]=ComputeBoundaryNodes(nodeweights);

    dx=1/NZy;
    rand('seed',1)
    for i=1:length(interiornodes)
        refnodes(:,interiornodes(i))=refnodes(:,interiornodes(i))+2*jitter*dx*(0.5-rand(2,1));
    end
end

%Compute the reference mesh
refmesh=ComputeMesh(refnodes, topo);

%Define number of Nodes and Zones
NN=size(refnodes,2);
NZ=size(topo,2);

if strcmp(Method,'Q2Q1d')
    % Compute the center dof locations
    refcennodes=zeros(2,NZ);
    for i=1:NZ
        refcennodes(:,i)=ComputeCentroid(refmesh(:,:,i));
    end
    
    % Compute the edge values
    edgelist=ComputeEdgeList(topo);
    NE=size(edgelist,2);
    % Compute locations of edge nodes
    refedgenodes=(refnodes(:,edgelist(1,:))+refnodes(:,edgelist(2,:)))/2;
    
    % Construct the edge DOF mapping
    mapping=ComputeElementEdgeMapping2(topo,edgelist);
    
    % Define sub-cell mapping
    subcellMapping=[1 5 9 8; 5 2 6 9; 9 6 3 7; 8 9 7 4];
end

% Construct the Q1b map using the Node, Edge, and Cell map
Quadmap=zeros(ndofpZ,NZ);

for i=1:NZ
    % Start with the nodes
    Quadmap(1:4,i)=topo(:,i);
    if strcmp(Method,'Q2Q1d')
        % Add the edges
        Quadmap(5:8,i)=NN+mapping(:,i);
        % Add the centers
        Quadmap(9,i)=NN+NE+i;
    end
end

% Compute nodal and zonal strides
NSTRIDE=NZx+1;
ZSTRIDE=NZx;

% Compute number of dofs in Domain
if strcmp(Method,'Q2Q1d')
    ndofpD=NN+NE+NZ;
else
    ndofpD=NN;
end

% Construct the boundary edge sets
[bdof nbnodes]=ComputeBoundaryDOFQuad(Quadmap,isEdgeBasis,ZSTRIDE);