function refmesh=ComputeMesh(refnodes, topo)

refmesh=zeros(2,4,size(topo,2));

for i=1:size(topo,2)
    for j=1:4
        refmesh(:,j,i)=refnodes(:,topo(j,i));
    end
end