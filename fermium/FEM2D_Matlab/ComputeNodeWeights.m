function nodeweights=ComputeNodeWeights(topo)

NN=max(max(topo));
NZ=size(topo,1);
nodeweights=zeros(NN,1);

for i=1:NZ
    for j=1:4
        nodeweights(topo(i,j))=nodeweights(topo(i,j))+1;
    end
end