% function massM=GetMassMatrix(nodes,pts2d,wgts2d)

% clear
clc

[points,weights] = GaussLobattoWeights(3);
wgts2d=zeros(length(weights)^2,1);
pts2d=zeros(length(points)^2,2);
for i=1:length(weights)
    for j=1:length(weights)
        wgts2d(j+(i-1)*length(weights))=weights(i)*weights(j);
        pts2d(j+(i-1)*length(points),1)=points(i);
        pts2d(j+(i-1)*length(points),2)=points(j);
    end
end
nodes(:,:)=pertmesh(1,:,:);

Jmat=zeros(size(pts2d,1),2,2);
detJ=zeros(size(pts2d,1),1);


for n=1:size(pts2d,1)
    Jmat(n,:,:)=JacobianMat(nodes,pts2d(n,:));
    J(:,:)=Jmat(n,:,:);
    detJ(n)=det(J);
end
Jmat
detJ

massM=zeros(4);

for i=1:4
    for j=1:4
        