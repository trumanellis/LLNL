%% Work and EOS   
% Calculate the new zone volumes from the new node coordinates
if ~isPresRefined
   for i=1:NZ
       NEWvolumeZ(i)=ComputeVolume(newmesh(:,:,i));
       NEWdensityZ(i)=massZ(i)/NEWvolumeZ(i);

       % Compute the new internal energy
       % F dot V internal energy update
       zonevel=0.5*(NEWvelocityQuad(:,Quadmap(:,i))+OLDvelocityQuad(:,Quadmap(:,i)));
       NEWenergyZ(i)=OLDenergyZ(i)-dt/massZ(i)*sum(dot(zonevel,cornerForce(:,:,i)));
   end
else
   for i=1:NZ
       % Zonal Quantities
       % Compute the new internal energy
       % F dot V internal energy update
       zonevel=0.5*(NEWvelocityQuad(:,Quadmap(:,i))+OLDvelocityQuad(:,Quadmap(:,i)));
       NEWenergyZ(i)=OLDenergyZ(i)-dt/massZ(i)*sum(dot(zonevel,cornerForce(:,:,i)));

       % Sub-Zonal Quantities
       zonenodes=allnodes(:,Quadmap(:,i));
       for k=1:4
           NEWsubvolume(k,i)=ComputeVolume(zonenodes(:,subcellMapping(:,k)));
       end
       NEWdensityQ1d(:,i)=cornerMass(:,i)./NEWsubvolume(:,i);
       NEWdensityZ(i)=mean(NEWdensityQ1d(:,i));
   end
end

% ===EOS PHASE ===

% Evaluate the EOS for an ideal monatomic gas to get the new pressure
if ~isPresRefined
   for i=1:NZ
       pressureZ(i)=(gamma-1)*NEWdensityZ(i)*NEWenergyZ(i);
       if pressureZ(i) < 0
           fprintf('WARNING -- Negative Pressure detected at cell %4.0f \n',i)
       end
   end
else
   for i=1:NZ
       % Sub-Zonal Quantities
       pressureQ1d(:,i)=(gamma-1)*NEWdensityQ1d(:,i)*NEWenergyZ(i);
       pressureZ(i)=mean(pressureQ1d(:,i));
       for s=1:4
           if pressureQ1d(s,i) < 0
               pressureQ1d(s,i)=0;
               fprintf('WARNING -- Negative Pressure detected at cell %4.0f \n',i)
           end
       end
   end
end

% Tally up the energies (This is actually the energy for the previous cycle)
internalenergy(cycle)=sum(massZ.*OLDenergyZ);
if ~isFullMassMatrixSolve
   kineticenergy(cycle)=0.5*(OLDvelocityQuad(1,:)*MMQuad*OLDvelocityQuad(1,:)'+OLDvelocityQuad(2,:)*MMQuad*OLDvelocityQuad(2,:)');
   sumMassN(cycle)=sum(massN(:));
else
   kineticenergy(cycle)=0.5*(OLDvelocityQuad(1,:)*MMX*OLDvelocityQuad(1,:)'+OLDvelocityQuad(2,:)*MMY*OLDvelocityQuad(2,:)');
   sumMassN(cycle)=sum(MMQuad(:));
end

% Increment time
t=t+dt;
if t>=tstop
   stop=true;
end
timedata(cycle)=t;

% Reset old state variables to new state variables
OLDvelocityQuad=NEWvelocityQuad;
OLDvolumeZ=NEWvolumeZ;
OLDdensityZ=NEWdensityZ;
OLDenergyZ=NEWenergyZ;
if strcmp(Method,'Q2Q1d')
   oldedgenodes=newedgenodes;
   oldcennodes=newcennodes;
   OLDdensityQ1d=NEWdensityQ1d;
   OLDsubvolume=NEWsubvolume;
end

oldmesh=newmesh;
oldnodes=newnodes;