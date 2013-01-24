% QuadQ0 Acoustic Wave

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=1;
% Set initial time increment
dtInit=5e-3;
dtMax=5e-3;

Method='Q2Q1d'
isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=false;
massQuadOrder=3;
stiffQuadOrder=3;

Qfrac=0;

NZx=9;
NZy=9;

xmin=0;
xmax=1;

ymin=0;
ymax=xmax;

% We will run the problem to a final time of 0.9
tstart = 0;
tstop = .4;
vsound = 1;

%Define a random "jittering" factor -- this determines the magnitude of
%mesh distortion
jitter=0.2;

switch Method
    case 'Q1Q0_hgOFF'
        Basis='Q1Basis'
        GradBasis='GradQ1Basis'
        % Number of dofs per Zone
        ndofpZ=4;
        isEdgeBasis=false;
        isCenterBasis=false;
        isVelRefined=false;
        isPresRefined=false;
        hgfrac=0;
    case 'Q1Q0_hgON'
        Basis='Q1Basis'
        GradBasis='GradQ1Basis'
        % Number of dofs per Zone
        ndofpZ=4;
        isEdgeBasis=false;
        isCenterBasis=false;
        isVelRefined=false;
        isPresRefined=true;
        hgfrac=1;
    case 'Q2Q1d'
        Basis='Q2Basis'
        GradBasis='GradQ2Basis'
        % Number of dofs per Zone
        ndofpZ=9;
        isEdgeBasis=true;
        isCenterBasis=true;
        isVelRefined=false;
        isPresRefined=true;
        hgfrac=0;
    case 'Q1Q0'
        Basis='Q1Basis'
        GradBasis='GradQ1Basis'
        % Number of dofs per Zone
        ndofpZ=4;
        isEdgeBasis=false;
        isCenterBasis=false;
        isVelRefined=false;
        isPresRefined=false;
        hgfrac=0;
end

if SaveFigures
    if ~isVelRefined && ~isPresRefined
        figurefile=['FiguresInitGausPres_',Basis,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    elseif isVelRefined
        figurefile=['FiguresInitGausPres_Refined',Basis,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    elseif isPresRefined
        figurefile=['FiguresInitGausPres_Q1dPres',Basis,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    end
    unix(['mkdir FigureFiles/',figurefile]);
end

% Specify the maximum number of time steps to use
maxcycle=1;

% These numbers are the coefficients for the artificial viscosity
qquad=2;
qlin=0.25;

%% Mesh Assembly
MeshAssembly

%% Variable Preallocation
Preallocate

%% Initialize Properties
% Set the EOS gamma law constant for an ideal gas
gamma=5/3;
% pressFun=inline('x(1).*(1-x(1))','x');
pressFun=inline('sin(pi*x(1))','x');
rhoInit=1;

% Initialize the Zonal state variables
OLDdensityZ(:)=rhoInit;
if strcmp(Method,'Q2Q1d')
    OLDdensityQ1d(:)=rhoInit;
end

for i=1:NZ
    cen=ComputeCentroid(refmesh(:,:,i));
    pressureZ(i)=pressFun(cen);
    if strcmp(Method,'Q2Q1d')
        pressureQ1d(1,i)=pressFun(LocalToGlobal(refmesh(:,:,i),[.25 .25]));
        pressureQ1d(2,i)=pressFun(LocalToGlobal(refmesh(:,:,i),[.75 .25]));
        pressureQ1d(3,i)=pressFun(LocalToGlobal(refmesh(:,:,i),[.75 .75]));
        pressureQ1d(4,i)=pressFun(LocalToGlobal(refmesh(:,:,i),[.25 .75]));
    end
end
OLDenergyZ=pressureZ/((gamma-1)*rhoInit);

for i=1:NZ
    OLDvolumeZ(i)=ComputeVolume(refmesh(:,:,i));
    massZ(i)=OLDdensityZ(i)*OLDvolumeZ(i);
    if strcmp(Method,'Q2Q1d')
        MMzone=GetMassMatrix(refmesh(:,:,i),'Q1Basis',2);
        cornerMass(:,i)=MMzone*OLDdensityQ1d(:,1);
    end
end

% Assemble Nodal Masses if Mass is not updated every cycle
compmesh=refmesh;
if ~isMassUpdateEveryCycle && ~isFullMassMatrixSolve
    NodalMassAssembly
end

%% === Generation Visualization Phase ===

% Step 2) Now we loop over time steps and solve the hydro equations at each
% Lagrangian time step

% === Lagrangian Phase ===

%% === Time Stepping Loop ===
for cycle=1:maxcycle
    if stop==true
        break;
    end
    
    % Compute stable time increment
    StableTimeStep
    
%% === Acceleration Phase ===
    % Compute Corner Forces
    CornerForces

    if isFullMassMatrixSolve
        %% Mass Matrix Assembly
        MassMatrixAssembly
        
        % Impose Boundary Conditions
        % Impose symmetric boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(2,index)=0;
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(2,index)=0;
        end

        % Solve Linear System
        accelerationQuad=[(MMX\forceN(1,:)')'; (MMY\forceN(2,:)')'];
        
    elseif isMassUpdateEveryCycle
        % Assemble Nodal Masses
        compmesh=oldmesh;
        NodalMassAssembly
        
        % Solve for acceleration
        for i=1:ndofpD
            accelerationQuad(:,i)=forceN(:,i)/massN(i);
        end
        
        accelerationTransport(1,:)=forceTransport(1,:)./massTransport;
        accelerationTransport(2,:)=forceTransport(2,:)./massTransport;
        
        % Impose Boundary Conditions
        % Impose symmetric boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            accelerationQuad(2,index)=0;
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            accelerationQuad(2,index)=0;
        end
    else
        % Solve for acceleration
        for i=1:ndofpD
            accelerationQuad(:,i)=forceN(:,i)/massN(i);
        end
        
        accelerationTransport(1,:)=forceTransport(1,:)./massTransport;
        accelerationTransport(2,:)=forceTransport(2,:)./massTransport;
        
        % Impose Boundary Conditions
        % Impose symmetric boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            accelerationQuad(2,index)=0;
            if index<=NN
                accelerationTransport(2,index)=0;
            end
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            accelerationQuad(2,index)=0;
            if index<=NN
                accelerationTransport(2,index)=0;
            end
        end
    end

    % Integrate accelerations to get velocities
    NEWvelocityQuad=OLDvelocityQuad+dt*accelerationQuad;
    NEWvelocityTransport=OLDvelocityTransport+dt*accelerationTransport;
    
    % Integrate velocities to get displacements and update the grid
   newnodes=oldnodes+dt*NEWvelocityTransport;
   newmesh=ComputeMesh(newnodes,topo);
   
%% === Work and EOS Phase ===   
   WorkandEOS
   
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
       OLDdensityQ1d=NEWdensityQ1d;
       OLDsubvolume=NEWsubvolume;
       OLDvelocityTransport=NEWvelocityTransport;
   end
   
   oldmesh=newmesh;
   oldnodes=newnodes;
   
   % === VISUALIZATION AND POST-PROCESSING PHASE ===
   
%% =========== PLOT RESULTS ===========
    if mod(cycle,plotcycle) == 0
        centerAcc=zeros(2,NZ);
        cen=zeros(2,NZ);
        for i=1:NZ
            zoneacc=accelerationQuad(:,Quadmap(:,i));
            cen(:,i)=ComputeCentroid(refmesh(:,:,i));
            centerAcc(:,i)=LocalInterpolate(zoneacc,[0.5 0.5],Basis);
        end
       figure(1)
       hold on
       if strcmp(Method,'Q2Q1d')
           plot(cen(1,:),accelerationQuad(1,Quadmap(9,:)),'mo')
       else
           plot(cen(1,:),centerAcc(1,:),'mo')
       end
        plot(refnodes(1,:),accelerationQuad(1,1:NN),'bo')
        plot(refnodes(1,:),accelerationTransport(1,1:NN),'go')
        xlabel('x-position')
        ylabel('Acceleration')
        
        figure(2)
       ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),pressureZ);
       hold on
       nplotpoints=3;
       h=PlotFEMVectorFieldQ2(refmesh,accelerationQuad/20,Quadmap,nplotpoints,Basis,0,'m','LineWidth',1);
       axis equal
       axis tight
    end  
end