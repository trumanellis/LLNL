%% Variable Preallocation

% Initialize Quad state variables
accelerationQuad=zeros(2,ndofpD);
OLDvelocityQuad=zeros(2,ndofpD);
NEWvelocityQuad=zeros(2,ndofpD);

% Initialize the Zonal state variables
OLDenergyZ=zeros(1,NZ);
NEWenergyZ=zeros(1,NZ);

OLDdensityZ=zeros(1,NZ);
NEWdensityZ=zeros(1,NZ);

pressureZ=zeros(1,NZ);

OLDvolumeZ=zeros(1,NZ);
NEWvolumeZ=zeros(1,NZ);

massZ=zeros(1,NZ);
cornerForce=zeros(2,ndofpZ,NZ);

if strcmp(Method,'Q2Q1d')
    % Preallocate Q1d thermodynamic variables
    OLDdensityQ1d=zeros(4,NZ);
    pressureQ1d=zeros(4,NZ);
    NEWsubvolume=zeros(4,NZ);
    NEWdensityQ1d=zeros(4,NZ);
    cornerMass=zeros(4,NZ);
    
    % Initialize locations of "old" edge and center dofs
    oldedgenodes=refedgenodes;
    newedgenodes=refedgenodes;
    oldcennodes=refcennodes;
    newcennodes=refcennodes;
    allnodes=[refnodes,refedgenodes,refcennodes];
end

len=zeros(1,NZ);
cen=zeros(2,NZ);

% Initialize the "old" mesh
oldnodes=refnodes;
newnodes=refnodes;
oldmesh=refmesh;
newmesh=refmesh;

% Initialize time history arrays for integrated quantities
timedata=zeros(maxcycle,1);
internalenergy=zeros(maxcycle,1);
kineticenergy=zeros(maxcycle,1);
sumMassN=zeros(maxcycle,1);
cycle=0;

% Initialize looping conditions
stop=false;
t=tstart;