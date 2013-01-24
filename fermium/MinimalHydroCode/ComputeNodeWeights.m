function nodeweights=ComputeNodeWeights(topo)

NN=max(max(topo));
NZ=size(topo,2);
nodeweights=zeros(NN,1);

for i=1:NZ
    for j=1:4
        nodeweights(topo(j,i))=nodeweights(topo(j,i))+1;
    end
end