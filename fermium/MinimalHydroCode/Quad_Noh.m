% QuadQ Noh

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=5;
% Set initial and max time increment
dtInit=1e-4;
dtMax=5e-4;

% Define the method of choice
Method='Q2Q1d'
isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=true;
isEdgeCenIndependent=false;
massQuadOrder=3;
stiffQuadOrder=3;

% Use this to turn artificial viscosity on/off
Qfrac=1;

% Use this to turn anti-hourglass forces on/off (Only works with Q1Q0)
hgfrac=0;

% Define the number of zones
NZx=10;
NZy=10;

% Define the mesh size
xmin=0;
xmax=1;

ymin=0;
ymax=xmax;

% Define the start and stop time
tstart = 0;
tstop = 0.6;

%Define a random "jittering" factor -- this determines the magnitude of
%mesh distortion
jitter=0.0;

% Specify the maximum number of time steps to use
maxcycle=10000;

% These numbers are the coefficients for the artificial viscosity
qquad=2;
qlin=0.25;

switch Method
    case 'Q1Q0'
        Basis='Q1Basis'
        GradBasis='GradQ1Basis'
        % Number of dofs per Zone
        ndofpZ=4;
        isEdgeBasis=false;
        isCenterBasis=false;
        isVelRefined=false;
        isPresRefined=false;
    case 'Q2Q1d'
        Basis='Q2Basis'
        GradBasis='GradQ2Basis'
        % Number of dofs per Zone
        ndofpZ=9;
        isEdgeBasis=true;
        isCenterBasis=true;
        isVelRefined=false;
        isPresRefined=true;
end

if SaveFigures
    figurefile=['FiguresNoh_',Method,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    unix(['mkdir FigureFiles/',figurefile]);
end

%% Mesh Assembly
MeshAssembly

%% Variable Preallocation
Preallocate

%% Initialize Properties
% Set the EOS gamma law constant for an ideal gas
gamma=5/3;

rhoInit=1;
pInit=0;
eInit=0;

% Set initial velocity to radially inward
for i=1:NN
    thetaN=atan2(refnodes(2,i),refnodes(1,i));
    OLDvelocityQuad(:,i)=-1*[cos(thetaN); sin(thetaN)];
end
if strcmp(Basis,'Q2Basis')
    for i=1:NZ
        e1=(refmesh(:,1,i)+refmesh(:,2,i))/2;
        thetaN=atan2(e1(2),e1(1));
        OLDvelocityQuad(:,Quadmap(5,i))=-1*[cos(thetaN); sin(thetaN)];
        e2=(refmesh(:,2,i)+refmesh(:,3,i))/2;
        thetaN=atan2(e2(2),e2(1));
        OLDvelocityQuad(:,Quadmap(6,i))=-1*[cos(thetaN); sin(thetaN)];
        e3=(refmesh(:,3,i)+refmesh(:,4,i))/2;
        thetaN=atan2(e1(2),e1(1));
        OLDvelocityQuad(:,Quadmap(7,i))=-1*[cos(thetaN); sin(thetaN)];
        e4=(refmesh(:,4,i)+refmesh(:,1,i))/2;
        thetaN=atan2(e4(2),e4(1));
        OLDvelocityQuad(:,Quadmap(8,i))=-1*[cos(thetaN); sin(thetaN)];
        cen=ComputeCentroid(refmesh(:,:,i));
        thetaN=atan2(cen(2),cen(1));
        OLDvelocityQuad(:,Quadmap(9,i))=-1*[cos(thetaN); sin(thetaN)];
    end
end

% Set the velocity at the origin to be zero
OLDvelocityQuad(:,1)=[0;0];

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
    if strcmp(Method,'Q2Q1d')
        zonenodes=allnodes(:,Quadmap(:,i));
        for k=1:4
           NEWsubvolume(k,i)=ComputeVolume(zonenodes(:,subcellMapping(:,k)));
        end
        cornerMass(:,i)=NEWsubvolume(:,i).*OLDdensityQ1d(:,i);
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
    end

%% Update Mesh
   UpdateMesh
    
   
%% === Work and EOS Phase ===   
   WorkandEOS
   
%% === VISUALIZATION AND POST-PROCESSING PHASE ===  
   %% =========== PLOT RESULTS ===========
   if mod(cycle,plotcycle) == 0 || stop
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  Total Energy: %20.18f, Nodal Mass: %20.18f \n',cycle,t,internalenergy(cycle)+kineticenergy(cycle),sumMassN(cycle));
       plotvarZ=NEWdensityZ;
       plotvarN=sqrt(NEWvelocityQuad(1,:).^2+NEWvelocityQuad(2,:).^2);
       
       lineoutZ=plotvarZ(1:ZSTRIDE);
       lineoutN=plotvarN(1:NSTRIDE);
       [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ);
       
       sfigure(1);
       if strcmp(Method,'Q2Q1d')
           subplot(2,2,1)
           hold on
           subplot(2,2,2)
           hold on
           for i=1:NZ
               zonenodes=allnodes(:,Quadmap(:,i));
                zoneDOF=plotvarN(Quadmap(:,i));
                subcellMapping=[1 5 9 8; 5 2 6 9; 9 6 3 7; 8 9 7 4];
                subplot(2,2,1)
                trisurf2(subcellMapping,zonenodes(1,:),zonenodes(2,:),zoneDOF');

                subplot(2,2,2)
                ZonalSurf(subcellMapping,zonenodes(1,:),zonenodes(2,:),zeros(9,1),NEWdensityQ1d(:,i)');
           end
           subplot(2,2,1)
           axis equal tight
           title('Velocity Magnitude','FontWeight','demi','FontSize',14)
           az=0;
           el=90;
           view(az,el);
           colorbar
           hold off
       else
           subplot(2,2,1)
           trisurf2(topo',newnodes(1,:),newnodes(2,:),plotvarN');
           axis equal tight
           title('Velocity Magnitude','FontWeight','demi','FontSize',14)
           az=0;
           el=90;
           view(az,el);
           colorbar

           subplot(2,2,2)
           ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),plotvarZ);
           axis equal tight
           title('Density','FontWeight','demi','FontSize',14)
           colorbar
       end
       
       subplot(2,2,2)
%        ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),plotvarZ);
       axis equal tight
       title('Density','FontWeight','demi','FontSize',14)
       colorbar
       hold off
       
       subplot(2,2,3)
       hold on
       plot(timedata(1:cycle),(internalenergy(1:cycle)+kineticenergy(1:cycle))/(internalenergy(1)+kineticenergy(1)),'k')
       plot(timedata(1:cycle),internalenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'r')
       plot(timedata(1:cycle),kineticenergy(1:cycle)/(internalenergy(1)+kineticenergy(1)),'b')
       xlabel('Time [sec]','FontWeight','demi','FontSize',14)
       ylabel('Normalized Energy','FontWeight','demi','FontSize',14)
       legend('Total Energy','Internal Energy','Kinetic Energy','Location','NorthOutside')
       axis([0 t 0 1.2])

       subplot(2,2,4)
       plot(X(:,1),lineoutZ)
       title('Density','FontWeight','demi','FontSize',14)
       
       subplot_title(['Noh: t = ',num2str(t,'%5.4f')],[.1 .9 .8 .05]);
       
       if SaveFigures
           % Save each plot frame for conversion to animated gif
           saveas(1,['FigureFiles/',figurefile,'/plot',datestr(clock, 'ddHHMMSS'),'.png']);
       else
           drawnow
       end
   end
end
