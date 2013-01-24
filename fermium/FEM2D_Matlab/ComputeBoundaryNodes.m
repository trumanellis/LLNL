function [bnodes, interiornodes]=ComputeBoundaryNodes(nodeweights)

bc=1;
ic=1;
for i=1:length(nodeweights)
    if nodeweights(i)<4
        bnodes(bc)=i;
        bc=bc+1;
    else
        interiornodes(ic)=i;
        ic=ic+1;
    end
end

if ic==1
    interiornodes=[];
end