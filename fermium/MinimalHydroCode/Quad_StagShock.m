% QuadQ StagShock

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=10;
% Set initial time increment
dtInit=5e-4;
dtMax=5e-4;

Method='Q2Q1d'
isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=false;
massQuadOrder=3;
stiffQuadOrder=3;

Qfrac=1;

NZx=20;
NZy=1;

xmin=0;
xmax=1;

ymin=0;
ymax=0.1*xmax;

% We will run the problem to a final time of 0.9
tstart = 0;
tstop = 0.9;

%Define a random "jittering" factor -- this determines the magnitude of
%mesh distortion
jitter=0.0;

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
    figurefile=['FiguresStagShock_',Method,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    unix(['mkdir FigureFiles/',figurefile]);
end

% Specify the maximum number of time steps to use
maxcycle=5000;

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

rhoInit=1;
pInit=0;
eInit=1e-6;

% Set initial velocity to be uniform
OLDvelocityQuad(1,:)=1;
OLDvelocityQuad(:,bdof(1+nbnodes(1):sum(nbnodes(1:2))))=0;

% Initialize the Zonal state variables
OLDenergyZ(:)=eInit;
OLDdensityZ(:)=rhoInit;
pressureZ(:)=pInit;

if strcmp(Method,'Q2Q1d')
    OLDdensityQ1d(:)=rhoInit;
    pressureQ1d(:)=pInit;
end

for i=1:NZ
    OLDvolumeZ(i)=ComputeVolume(refmesh(:,:,i));
    massZ(i)=OLDdensityZ(i)*OLDvolumeZ(i);
end
for i=1:NZ
    OLDvolumeZ(i)=ComputeVolume(refmesh(:,:,i));
    massZ(i)=OLDdensityZ(i)*OLDvolumeZ(i);
    if strcmp(Method,'Q2Q1d')
        MMzone=GetMassMatrix(refmesh(:,:,i),'Q1Basis',2);
        cornerMass(:,i)=MMzone*OLDdensityQ1d(:,i);
    end
end

% Assemble Nodal Masses if Mass is not updated every cycle
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
        % Impose zero boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
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
        NodalMassAssembly
        
        % Solve for acceleration
        for i=1:ndofpD
            accelerationQuad(:,i)=forceN(:,i)/massN(i);
        end
        
        % Impose Boundary Conditions
        % Impose symmetric boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            accelerationQuad(2,index)=0;
        end
        % Impose zero boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
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
        
        % Impose Boundary Conditions
        % Impose zero boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            accelerationQuad(2,index)=0;
        end
        % Impose zero boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
        end
        % Impose symmetric boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            accelerationQuad(2,index)=0;
        end
    end

    % Integrate accelerations to get velocities
    NEWvelocityQuad=OLDvelocityQuad+dt*accelerationQuad;

    % Integrate velocities to get displacements and update the grid
   newnodes=oldnodes+dt*NEWvelocityQuad(:,1:NN);
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
   end
   
   oldmesh=newmesh;
   oldnodes=newnodes;
   
%% =========== PLOT RESULTS ===========
   if mod(cycle,plotcycle) == 0
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  Total Energy: %20.18f, Nodal Mass: %20.18f \n',cycle,t,internalenergy(cycle)+kineticenergy(cycle),sumMassN(cycle));
       plotvarZ=pressureZ;
       plotvarN=NEWvelocityQuad(1,1:NN);
     
       lineoutZ2=pressureZ(1:ZSTRIDE);
       lineoutZ1=NEWdensityZ(1:ZSTRIDE);
       lineoutN=plotvarN(1:NSTRIDE);
       [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ);
       
       sfigure(1);
       subplot(3,2,1:2)
       trisurf2(topo',newnodes(1,:),newnodes(2,:),NEWvelocityQuad(1,1:NN)');
       axis([xmin xmax ymin ymax])
       axis equal
       title('x-Velocity','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
       subplot(3,2,3)
       ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),NEWdensityZ);
       axis tight
       title('Density','FontWeight','demi','FontSize',14)
       colorbar
       
       subplot(3,2,5)
       hold on
       plot(timedata(1:cycle),(internalenergy(1:cycle)+kineticenergy(1:cycle))/(internalenergy(1)+kineticenergy(1)),'k')
       plot(timedata(1:cycle),internalenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'r')
       plot(timedata(1:cycle),kineticenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'b')
       xlabel('Time [sec]','FontWeight','demi','FontSize',14)
       ylabel('Normalized Energy','FontWeight','demi','FontSize',14)
       legend('Total Energy','Internal Energy','Kinetic Energy','Location','NorthOutside')
       axis([0 t 0 1.2])

       subplot(3,2,[4 6])
       [AX,H1,H2]=plotyy(X(:,1),lineoutZ1,X(:,1),lineoutZ2);
%        set(AX(2),'YLim',[0 1.4]) 
%        set(AX(1),'YLim',[0 5])
%        set(AX(2),'YTick',linspace(0,1.4,5))
%        set(AX(2),'YTickLabel',num2str(linspace(0, 1.4,5)'))
%        set(AX(1),'YTick',linspace(0,5,5))
%        set(AX(1),'YTickLabel',num2str(linspace(0, 5,5)'))
%        set(AX(2),'YGrid','on')
       set(get(AX(2),'Ylabel'),'String','Pressure','FontWeight','demi','FontSize',10) 
       set(get(AX(1),'Ylabel'),'String','Density','FontWeight','demi','FontSize',10)       
       
       subplot_title(['Stagnation Shock: t = ',num2str(t)],[.1 .9 .8 .05]);
       
       if SaveFigures
           % Save each plot frame for conversion to animated gif
           saveas(1,['FigureFiles/',figurefile,'/plot',datestr(clock, 'ddHHMMSS'),'.png']);
       else
           drawnow
       end
   end
   
end