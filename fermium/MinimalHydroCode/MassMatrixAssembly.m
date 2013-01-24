%% Full Mass Matrix Assembly
MMQuad=zeros(ndofpD);

for N=1:NZ
    if ~isPresRefined
        localm=OLDdensityZ(N)*GetMassMatrix(oldmesh(:,:,N),Basis,massQuadOrder);
    else
        localm=GetMassMatrix(oldmesh(:,:,N),Basis,massQuadOrder,OLDdensityQ1d(:,N));
    end
    MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+localm;
end
MMX=MMQuad;
MMY=MMQuad;