function stiffM=GetStiffnessMatrix(nodes,pts2d,wgts2d)
% clc
% 
% [points,weights] = GaussLobattoWeights(3);
% wgts2d=zeros(length(weights)^2,1);
% pts2d=zeros(length(points)^2,2);
% for i=1:length(weights)
%     for j=1:length(weights)
%         wgts2d(j+(i-1)*length(weights))=weights(i)*weights(j);
%         pts2d(j+(i-1)*length(points),1)=points(i);
%         pts2d(j+(i-1)*length(points),2)=points(j);
%     end
% end
% nodes(:,:)=pertmesh(1,:,:);

Jmat=zeros(2,2);
detJ=zeros(1);
invJ=zeros(2,2);

stiffM=zeros(4);
for n=1:size(pts2d,1)
    Jmat(:,:)=JacobianMat(nodes,pts2d(n,:));
    detJ=det(Jmat);
    invJ=inv(Jmat);
    
    [dwdx dwdy]=BasisDeriv(pts2d(n,1),pts2d(n,2));
    stiffM=stiffM+wgts2d(n)*[[dwdx dwdy]*invJ]*[[dwdx dwdy]*invJ]'*detJ;
end

