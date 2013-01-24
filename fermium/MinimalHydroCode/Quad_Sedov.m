% QuadQ0 Acoustic Wave

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=1;
% Set initial time increment
dtInit=1e-5;
dtMax=1e-4;

Method='Q2Q1d'
isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=false;
massQuadOrder=3;
stiffQuadOrder=3;

Qfrac=1;

NZx=10;
NZy=10;

xmin=0;
xmax=1;

ymin=0;
ymax=xmax;

% We will run the problem to a final time of 0.9
tstart = 0;
tstop = 0.1;

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
        isPresRefined=false;
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
    figurefile=['FiguresSedov_',Method,datestr(clock, 'yyyy_mm_dd_HH_MM')];
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
eInit=1e-6;
eBigInit=1e2;
pInit=(gamma-1)*rhoInit*eInit;
pBigInit=(gamma-1)*rhoInit*eBigInit;

% Initialize the Zonal state variables
OLDenergyZ(:)=eInit;
OLDenergyZ(1)=eBigInit;
OLDdensityZ(:)=rhoInit;
pressureZ(:)=pInit;
pressureZ(1)=pBigInit;

if strcmp(Method,'Q2Q1d')
    OLDdensityQ1d(:)=rhoInit;
    pressureQ1d(:)=pInit;
    pressureQ1d(:,1)=pBigInit;
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
        % Impose symmetric boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            forceN(1,index)=0;
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
        % Impose symmetric boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            accelerationQuad(1,index)=0;
        end
    else
        % Solve for acceleration
        for i=1:ndofpD
            accelerationQuad(:,i)=forceN(:,i)/massN(i);
        end
        
%         accelerationTransport(1,:)=forceTransport(1,:)./massTransport;
%         accelerationTransport(2,:)=forceTransport(2,:)./massTransport;
        
        % Impose Boundary Conditions
        % Impose symmetric boundary at bottom wall
        for i=1:sum(nbnodes(1))
            index=bdof(i);
            accelerationQuad(2,index)=0;
            if index<=NN
                accelerationTransport(2,index)=0;
            end
        end
        % Impose symmetric boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            accelerationQuad(1,index)=0;
            if index<=NN
                accelerationTransport(1,index)=0;
            end
        end
    end

    % Integrate accelerations to get velocities
    NEWvelocityQuad=OLDvelocityQuad+dt*accelerationQuad;
    velocityTransport=NEWvelocityQuad(:,1:NN);
    for i=1:NZ
        velocityTransport(:,Quadmap(1,i))=velocityTransport(:,Quadmap(1,i))+0.5*sum(NEWvelocityQuad(:,Quadmap([8 5],i)),2)+.25*NEWvelocityQuad(:,Quadmap([9],i));
        velocityTransport(:,Quadmap(2,i))=velocityTransport(:,Quadmap(2,i))+0.5*sum(NEWvelocityQuad(:,Quadmap([5 6],i)),2)+.25*NEWvelocityQuad(:,Quadmap([9],i));
        velocityTransport(:,Quadmap(3,i))=velocityTransport(:,Quadmap(3,i))+0.5*sum(NEWvelocityQuad(:,Quadmap([6 7],i)),2)+.25*NEWvelocityQuad(:,Quadmap([9],i));
        velocityTransport(:,Quadmap(4,i))=velocityTransport(:,Quadmap(4,i))+0.5*sum(NEWvelocityQuad(:,Quadmap([7 8],i)),2)+.25*NEWvelocityQuad(:,Quadmap([9],i));
    end
    velocityTransport=velocityTransport/4;
    
    % Integrate velocities to get displacements and update the grid
    newnodes=oldnodes+dt*velocityTransport;
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
   end
   
   oldmesh=newmesh;
   oldnodes=newnodes;
   
%% =========== PLOT RESULTS ===========
   if mod(cycle,plotcycle) == 0 ||stop
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  Total Energy: %20.18f, Nodal Mass: %20.18f \n',cycle,t,internalenergy(cycle)+kineticenergy(cycle),sumMassN(cycle));
       plotvarZ=NEWdensityZ;
       plotvarN=sqrt(NEWvelocityQuad(1,1:NN).^2+NEWvelocityQuad(2,1:NN).^2);
       
       lineoutZ=plotvarZ(1:ZSTRIDE);
       lineoutN=plotvarN(1:NSTRIDE);
       [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ);
       
       sfigure(1);
       subplot(2,2,1)
       trisurf2(topo',newnodes(1,:),newnodes(2,:),plotvarN');
       axis equal tight
       title('Velocity Magnitude','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
       subplot(2,2,2)
       plot(newnodes(1,1:NSTRIDE),lineoutN);
       title('Velocity Magnitude','FontWeight','demi','FontSize',14)
       axis([xmin xmax 0 10])
       
       subplot(2,2,3)
       ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),plotvarZ);
       axis equal tight
       title('Density','FontWeight','demi','FontSize',14)
       colorbar

       subplot(2,2,4)
       hold on
       plot(timedata(1:cycle),(internalenergy(1:cycle)+kineticenergy(1:cycle))/(internalenergy(1)+kineticenergy(1)),'k')
       plot(timedata(1:cycle),internalenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'r')
       plot(timedata(1:cycle),kineticenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'b')
       xlabel('Time [sec]','FontWeight','demi','FontSize',14)
       ylabel('Normalized Energy','FontWeight','demi','FontSize',14)
       legend('Total Energy','Internal Energy','Kinetic Energy','Location','NorthOutside')
       axis([0 t 0 1.2])
       
       subplot_title(['Sedov: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
       
       if SaveFigures
           % Save each plot frame for conversion to animated gif
           saveas(1,['FigureFiles/',figurefile,'/plot',datestr(clock, 'ddHHMMSS'),'.png']);
       else
           drawnow
       end
   end
   
end

% quit