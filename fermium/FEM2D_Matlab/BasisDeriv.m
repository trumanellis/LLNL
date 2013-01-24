function [dwdx dwdy]=BasisDeriv(x,y)

dwdx=[
    -1+y;
    1-y;
    y;
    -y];

dwdy=[
    -1+x;
    -x;
    x;
    1-x];