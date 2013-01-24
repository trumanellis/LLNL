function [w]=Q2Basis(x,y)

w=zeros(9,length(x));

w(1,:)=(1+(-2+2*(-0.5+x)).*x).*(1+(-2+2*(-0.5+y)).*y);
w(2,:)=2*(-0.5+x).*x.*(1+(-2+2*(-0.5+y)).*y);
w(3,:)=4*(-0.5+x).*x.*(-0.5+y).*y;
w(4,:)=2*(1+(-2+2*(-0.5+x)).*x).*(-0.5+y).*y;
w(5,:)=(2-4*(-0.5+x)).*x.*(1+(-2+2*(-0.5+y)).*y);
w(6,:)=2*(-0.5+x).*x.*(2-4*(-0.5+y)).*y;
w(7,:)=2*(2-4*(-0.5+x)).*x.*(-0.5+y).*y;
w(8,:)=(1+(-2+2*(-0.5+x)).*x).*(2-4*(-0.5+y)).*y;
w(9,:)=(2-4*(-0.5+x)).*x.*(2-4*(-0.5+y)).*y;