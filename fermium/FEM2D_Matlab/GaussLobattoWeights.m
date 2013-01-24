function [points,weights] = GaussLobattoWeights(order)

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