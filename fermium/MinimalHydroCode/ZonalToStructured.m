function [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ)

X=zeros(NZx,NZy);
Y=zeros(NZx,NZy);
plotvarS=zeros(NZx,NZy);
for j=1:NZy
    for i=1:NZx
        index=i+(j-1)*NZx;
        nodes(:,:)=newmesh(:,:,index);
        temp=ComputeCentroid(nodes);
        X(i,j)=temp(1);
        Y(i,j)=temp(2);
        plotvarS(i,j)=plotvarZ(index);
    end
end