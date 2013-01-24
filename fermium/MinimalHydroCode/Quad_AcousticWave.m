% QuadQ0 Acoustic Wave

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=1;
% Set initial time increment
dtInit=1e-2;
dtMax=1e-2;

Method='Q2Q1d'
isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=false;
massQuadOrder=3;
stiffQuadOrder=3;

% Use this to turn artificial viscosity on/off
Qfrac=0;

% Use this to turn anti-hourglass forces on/off (Only works with Q1Q0)
hgfrac=0;

% Define the number of zones
NZx=20;
NZy=10;

% Define the mesh size
xmin=-1;
xmax=1;

ymin=0;
ymax=xmax;

% Define the start and stop time
tstart = 0;
tstop = 1;

%Define a random "jittering" factor -- this determines the magnitude of
%mesh distortion
jitter=0.0;

% Specify the maximum number of time steps to use
maxcycle=5000;

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
    figurefile=['FiguresAcousticWave_',Method,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    unix(['mkdir FigureFiles/',figurefile]);
end

%% Mesh Assembly
MeshAssembly

%% Variable Preallocation
Preallocate

%% Initialize Properties
% Set the EOS gamma law constant for an ideal gas
gamma=5/3;
vsound=1;

rhoInit=1;
pInit=vsound^2*rhoInit/gamma;
eInit=pInit/((gamma-1)*rhoInit);

amp=1/100;
omega=3*2*pi;
OLDvelocityQuad(1,floor(NSTRIDE/2)+1)=amp*cos(omega*tstart);

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
        % Impose zero boundary at bottom wall
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
        % Impose zero boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
        end
        % Impose symmetric boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            MMX(index,:)=0;
            MMX(:,index)=0;
            MMX(index,index)=1;
            MMY(index,:)=0;
            MMY(:,index)=0;
            MMY(index,index)=1;
            forceN(:,index)=[0;0];
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
        % Impose zero boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
        end
        % Impose zero boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
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
        % Impose zero boundary at right wall
        for i=1+nbnodes(1):sum(nbnodes(1:2))
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
        end
        % Impose zero boundary at top wall
        for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
        end
        % Impose zero boundary at left wall
        for i=1+sum(nbnodes(1:3)):sum(nbnodes)
            index=bdof(i);
            accelerationQuad(:,index)=[0;0];
        end
    end

%% Update Mesh
   UpdateMesh
   % ----- Apply Source term to velocities -----
    % This is where we "wiggle" the center node by specifying its velocity
    % for all time
    NEWvelocityQuad(1,floor(NSTRIDE/2)+1)=amp*cos(omega*t);
   
%% === Work and EOS Phase ===   
   WorkandEOS
   
%% === VISUALIZATION AND POST-PROCESSING PHASE ===
   
    %% =========== PLOT RESULTS ===========
   if mod(cycle,plotcycle) == 0 || stop
       fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  Total Energy: %20.18f, Nodal Mass: %20.18f \n',cycle,t,internalenergy(cycle)+kineticenergy(cycle),sumMassN(cycle));
       plotvarZ=pressureZ;
       plotvarN=NEWvelocityQuad(1,1:NN);
     
       lineoutZ=plotvarZ(1:ZSTRIDE);
       lineoutN=plotvarN(1:NSTRIDE);
       [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ);
       
       sfigure(1);
       subplot(2,2,1)
       trisurf2(topo',newnodes(1,:),newnodes(2,:),plotvarN');
       axis([xmin xmax ymin ymax -amp amp -1e-4 1e-4])
       axis equal tight
       title('x-Velocity','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
       subplot(2,2,3)
%        ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),plotvarZ,[.5995 .6005]);
       contourf(X,Y,plotvarS,[.59,linspace(.599,.601,20),.61],'LineStyle','none','LineWidth',0.5);
       axis equal tight
       caxis([.599 .601])
       title('Pressure','FontWeight','demi','FontSize',14)
       colorbar

       subplot(2,2,[2 4])
       [AX,H1,H2]=plotyy(X(:,1),lineoutZ,newnodes(1,1:NSTRIDE),lineoutN);
       set(AX(1),'YLim',[.598 .602]) 
       set(AX(2),'YLim',[-amp amp])
       set(AX(1),'YTick',linspace(.598,.602,5))
       set(AX(1),'YTickLabel',num2str(linspace(.598, .602,5)'))
       set(AX(2),'YTick',linspace(-amp,amp,5))
       set(AX(2),'YTickLabel',num2str(linspace(-amp, amp,5)'))
       set(AX(1),'YGrid','on')
       set(get(AX(1),'Ylabel'),'String','Pressure','FontWeight','demi','FontSize',10) 
       set(get(AX(2),'Ylabel'),'String','x-Velocity','FontWeight','demi','FontSize',10)       
       
       subplot_title(['Acoustic Wave: t = ',num2str(t,'%4.3f')],[.1 .85 .4 .05]);
       
       if SaveFigures
           % Save each plot frame for conversion to animated gif
           saveas(1,['FigureFiles/',figurefile,'/plot',datestr(clock, 'ddHHMMSS'),'.png']);
       else
           drawnow
       end
   end
   
end