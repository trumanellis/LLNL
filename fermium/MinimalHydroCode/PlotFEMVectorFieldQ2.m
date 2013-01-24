function h=PlotFEMVectorFieldQ2(mesh,projvector,mapping,npts,Basis,varargin)

temppts=linspace(0,1,npts);
localpts=zeros(2,size(temppts,2)^2);
for i=1:length(temppts)
    for j=1:length(temppts)
        localpts(1,j+(i-1)*length(temppts))=temppts(i);
        localpts(2,j+(i-1)*length(temppts))=temppts(j);
    end
end

pts=[];
vals=[];

globalpts=zeros(2,size(localpts,2));
globalvals=zeros(2,size(localpts,2));
for N=1:size(mesh,3)
    for i=1:size(localpts,2)
        globalpts(:,i)=LocalToGlobal(mesh(:,:,N),localpts(:,i));
    end
    vecproj=projvector(:,mapping(:,N));
    for i=1:size(localpts,2) 
        globalvals(:,i)=LocalInterpolateVelocityQ2(vecproj,localpts(:,i),Basis);
    end
    pts=[pts,globalpts];
    vals=[vals,globalvals];
end

h=quiver(pts(1,:),pts(2,:),vals(1,:),vals(2,:),varargin{:});