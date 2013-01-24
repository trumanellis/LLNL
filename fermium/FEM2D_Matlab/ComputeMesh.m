function refmesh=ComputeMesh(refnodes, topo)

refmesh=zeros(size(topo,1),4,2);

for i=1:size(topo,1)
    for j=1:4
        refmesh(i,j,:)=refnodes(topo(i,j),:);
    end
end