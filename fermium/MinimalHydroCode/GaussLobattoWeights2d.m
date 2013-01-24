function [pts2d,wgts2d] = GaussLobattoWeights2d(order)

switch order
    case 1
        points=[0,0];
        weights=[2];
    case 2
        points = [-1.0, 1.0];
        weights = [1.0, 1.0];
    case 3
        points = [-1.0, 0.0, 1.0];
        weights = [0.333333, 1.333333, 0.333333];
    case 4
        points = [-1.0, -0.447214, 0.447214, 1.0];
        weights = [0.166667, 0.833333, 0.833333, 0.166667];
    case 5
        points = [-1.0, -0.654654, 0.0, 0.654654, 1.0];
        weights = [0.1, 0.544444, 0.711111, 0.544444, 0.1];
end

%Reshift the interval to 0,1
weights=weights*.5;
points=0.5+0.5*points;

wgts2d=zeros(length(weights)^2,1);
pts2d=zeros(size(points,1)^2,2);
for i=1:length(weights)
    for j=1:length(weights)
        wgts2d(j+(i-1)*length(weights))=weights(i)*weights(j);
        pts2d(j+(i-1)*length(points),1)=points(i);
        pts2d(j+(i-1)*length(points),2)=points(j);
    end
end