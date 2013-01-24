%% Nodal Mass Assembly
% Nodal Masses
massN=zeros(ndofpD,1);
% Diagonal mass matrix for kinetic energy calculation
MMQuad=zeros(ndofpD);
cornerMassQ2=zeros(9,9,NZ);

for N=1:NZ
    if ~isPresRefined
        localm=OLDdensityZ(N)*GetMassMatrix(oldmesh(:,:,N),Basis,massQuadOrder);
        massN(Quadmap(:,N))=massN(Quadmap(:,N))+diag(localm,0);
        MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+localm;
    else
        cornerMassQ2(:,:,N)=GetMassMatrix(oldmesh(:,:,N),Basis,massQuadOrder,OLDdensityQ1d(:,N));
        massN(Quadmap(:,N))=massN(Quadmap(:,N))+diag(cornerMassQ2(:,:,N),0);
        MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+cornerMassQ2(:,:,N);
     end
end