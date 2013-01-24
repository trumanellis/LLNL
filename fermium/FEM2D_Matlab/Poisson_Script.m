%Using mixed finite elements to compute accelerations

%Step 1) 
%Initialize the mesh and all of the necessary topological data structures

clear all;
close all;
clc;

% nelems=2;
NZx=20;
NZy=20;

%Define the mesh extents
xmin=0;
xmax=1;
ymin=0;
ymax=1;

%Compute the reference mesh nodes
refnodes=ComputeReferenceMeshNodes(NZx,NZy,[xmin,xmax],[ymin,ymax]);

%Compute the mesh connectivity
topo=ComputeMeshTopology(NZx,NZy);

%Compute the reference mesh
refmesh=ComputeMesh(refnodes, topo);

%Define number of Nodes and Zones
NN=length(refnodes);
NZ=length(topo);

%Construct the node weights: i.e. the numer of zones which share each node
nodeweights=ComputeNodeWeights(topo);

%Compute the boundary and interior nodes
[bnodes, interiornodes]=ComputeBoundaryNodes(nodeweights);

%Compute nodal and zonal strides
NSTRIDE=NZx+1;
ZSTRIDE=NZx;

%Perturb the interior nodes to generate a non-orthogonal mesh
pertnodes=refnodes;

%Define a random "jittering" factor -- this determines the magnitude of
%mesh distortion
jitter=0.2;

dx=1/NZy;

for i=1:length(interiornodes)
    pertnodes(interiornodes(i),:)=pertnodes(interiornodes(i),:)+2*jitter*dx*(0.5-rand(1,2));
end

%Now define the perturbed mesh
pertmesh=ComputeMesh(pertnodes,topo);

figure(1)
hold on
for i=1:NZ
    plot([pertmesh(i,1,1), pertmesh(i,2,1), pertmesh(i,3,1), pertmesh(i,4,1), pertmesh(i,1,1)],...
        [pertmesh(i,1,2), pertmesh(i,2,2), pertmesh(i,3,2), pertmesh(i,4,2), pertmesh(i,1,2)],'k')
end
axis equal


%STEP 2



%STEP 3 Mass and Stiffness Matrix Calculation and Assembly

%Define the order of the element quadrature rule
quadorder=3;

[points,weights] = GaussLobattoWeights(quadorder);
wgts2d=zeros(length(weights)^2,1);
pts2d=zeros(length(points)^2,2);
for i=1:length(weights)
    for j=1:length(weights)
        wgts2d(j+(i-1)*length(weights))=weights(i)*weights(j);
        pts2d(j+(i-1)*length(points),1)=points(i);
        pts2d(j+(i-1)*length(points),2)=points(j);
    end
end

%Mass Matrix Assembly
% MM=zeros(NN);

%Stiffness Matrix Assembly
SM=zeros(NN);
for n=1:NZ
    nodes(:,:)=pertmesh(n,:,:);
    locals=GetStiffnessMatrix(nodes,pts2d,wgts2d);
    for i=1:4
        for j=1:4
            SM(topo(n,i),topo(n,j)) = SM(topo(n,i),topo(n,j)) + locals(i,j);
        end
    end
end

g=zeros(NN,1);
for N=1:NZ
    localval=zeros(4,1);
    for i=1:length(pts2d)
        localval=localval+wgts2d(i)*myForcingFunction(pts2d(i,1),pts2d(i,2))*basisFunctions(pts2d(i,1),pts2d(i,2));
    end
    for i=1:4
        g(topo(N,i))=g(topo(N,i))+localval(i);
    end
end

for n=1:length(bnodes)
    g(bnodes(n))=0;
    SM(:,bnodes(n))=0;
    SM(bnodes(n),:)=0;
    SM(bnodes(n),bnodes(n))=1;
end

u=SM\g

ccsol=zeros(NZ,1);
for N=1:NZ
    for i=1:4
    localsol(i)=u(topo(N,i));
    end
    ccsol(N)=0;
    w=basisFunctions(.5,.5);
    for i=1:4
        ccsol(N)=ccsol(N)+localsol(i)*w(i);
    end
end

figure(2)
trisurf(topo,pertnodes(:,1),pertnodes(:,2),u)
    