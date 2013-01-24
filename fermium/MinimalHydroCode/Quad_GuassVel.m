% QuadQ0 Acoustic Wave

clear all
close all
clc
format compact
SaveFigures=false;
plotcycle=10;
% Set initial time increment
dtInit=1e-2;
dtMax=1e-2;

Basis='Q1bBasis'
GradBasis='GradQ1bBasis'
% Number of dofs per Zone
ndofpZ=8;
isEdgeBasis=true;
isCenterBasis=false;
isRefinedQuad=false;

isFullMassMatrixSolve=false;
isMassUpdateEveryCycle=false;
massQuadOrder=4;
stiffQuadOrder=2;

Qfrac=0;
hgfrac=0;

NZx=10;
NZy=10;

xmin=0;
xmax=1;

ymin=0;
ymax=xmax;

% We will run the problem to a final time of 0.9
tstart = 0;
tstop = 1;
vsound = 1;

if SaveFigures
    if ~isRefinedQuad
        figurefile=['FiguresAcousticWave_',Basis,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    elseif isRefinedQuad
        figurefile=['FiguresAcousticWave_Refined',Basis,datestr(clock, 'yyyy_mm_dd_HH_MM')];
    end
    unix(['mkdir FigureFiles/',figurefile]);
end

% Specify the maximum number of time steps to use
maxcycle=5000;

% These numbers are the coefficients for the artificial viscosity
qquad=2;
qlin=0.25;

%Compute the mesh nodes
refnodes=ComputeReferenceMeshNodes(NZx,NZy,[xmin,xmax],[ymin,ymax]);

%Compute the mesh connectivity
topo=ComputeMeshTopology(NZx,NZy);

%Compute the reference mesh
refmesh=ComputeMesh(refnodes, topo);

%Define number of Nodes and Zones
NN=size(refnodes,2);
NZ=size(topo,2);

if isEdgeBasis==true
    edgelist=ComputeEdgeList(topo);
    NE=size(edgelist,2);

    % Construct the edge DOF mapping
    mapping=ComputeElementEdgeMapping2(topo,edgelist);  % This takes very long computationally
end

% Construct the Q1b map using the Node, Edge, and Cell map
Quadmap=zeros(ndofpZ,NZ);

for i=1:NZ
    % Start with the nodes
    Quadmap(1:4,i)=topo(:,i);
    if isEdgeBasis==true
        % Add the edges
        Quadmap(5:8,i)=NN+mapping(:,i);
    end
    if isCenterBasis==true
        % Add the centers
        Quadmap(9,i)=NN+NE+i;
    end
end

% Compute nodal and zonal strides
NSTRIDE=NZx+1;
ZSTRIDE=NZx;

% Compute number of dofs in Domain
if isEdgeBasis == true
    if isCenterBasis == true
        ndofpD=NN+NE+NZ;
    else
        ndofpD=NN+NE;
    end
elseif isCenterBasis ==true
    ndofpD=NN+NZ;
else
    ndofpD=NN;
end

% Construct the boundary edge sets
[bdof nbnodes]=ComputeBoundaryDOFQuad(Quadmap,ndofpZ,ZSTRIDE);

%% Initialize Properties
% Set the EOS gamma law constant for an ideal gas
gamma=5/3;

rhoInit=1;
pInit=vsound^2*rhoInit/gamma;
eInit=pInit/((gamma-1)*rhoInit);

% Initialize Quad state variables
accelerationQuad=zeros(2,ndofpD);
OLDvelocityQuad=zeros(2,ndofpD);
NEWvelocityQuad=zeros(2,ndofpD);

pow=25;
coeff=0.25;
% for i=1:NN
% %     OLDvelocityQuad(:,i)=coeff*2*pow*refnodes(:,i)*exp(-pow*(refnodes(1,i)^2+refnodes(2,i)^2));
%     thetaN=atan2(refnodes(2,i),refnodes(1,i));
%     OLDvelocityQuad(:,i)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(refnodes(1,i)^2+refnodes(2,i)^2));
% end
% OLDvelocityQuad(:,1)=[0;0];
if (strcmp(Basis,'Q2Basis') | isRefinedQuad)
%     OLDvelocityQuad(1,:)=1;
    for i=1:NZ
        thetaN=atan2(refnodes(2,i),refnodes(1,i));
        OLDvelocityQuad(:,Quadmap(1:4,i))=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(refnodes(1,Quadmap(1:4,i)).^2+refnodes(2,Quadmap(1:4,i)).^2));
        e1=mean(refmesh(:,1:2,i));
        thetaN=atan2(e1(2),e1(1));
        OLDvelocityQuad(:,5)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(e1(1).^2+e1(2)^2));
        e2=mean(refmesh(:,2:3,i));
        thetaN=atan2(e2(2),e2(1));
        OLDvelocityQuad(:,6)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(e2(1).^2+e2(2)^2));
        e3=mean(refmesh(:,3:4,i));
        thetaN=atan2(e3(2),e3(1));
        OLDvelocityQuad(:,7)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(e3(1).^2+e3(2)^2));
        e4=mean(refmesh(:,[4 1],i));
        thetaN=atan2(e4(2),e4(1));
        OLDvelocityQuad(:,8)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(e4(1).^2+e4(2)^2));
        cen=mean(refmesh(:,1:4,i));
        thetaN=atan2(cen(2),cen(1));
        OLDvelocityQuad(:,9)=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(cen(1).^2+cen(2)^2));
    end
else
%     OLDvelocityQuad(1,1:NN)=1;
    for i=1:NZ
        thetaN=atan2(refnodes(2,i),refnodes(1,i));
        OLDvelocityQuad(:,Quadmap(1:4,i))=[cos(thetaN); sin(thetaN)]*coeff*exp(-pow*(refnodes(1,Quadmap(1:4,i)).^2+refnodes(2,Quadmap(1:4,i)).^2));
    end
end
plotvarN=sqrt(OLDvelocityQuad(1,1:NN).^2+OLDvelocityQuad(2,1:NN).^2);

% Initialize the Zonal state variables
OLDenergyZ=eInit*ones(1,NZ);
NEWenergyZ=zeros(1,NZ);

OLDdensityZ=rhoInit*ones(1,NZ);
NEWdensityZ=zeros(1,NZ);

pressureZ=pInit*ones(1,NZ);

OLDvolumeZ=zeros(1,NZ);
NEWvolumeZ=zeros(1,NZ);

cornerMass=zeros(4,NZ);
massZ=zeros(1,NZ);
cornerForce=zeros(2,ndofpZ,NZ);

len=zeros(1,NZ);
cen=zeros(2,NZ);

for i=1:NZ
    OLDvolumeZ(i)=ComputeVolume(refmesh(:,:,i));
    massZ(i)=OLDdensityZ(i)*OLDvolumeZ(i);
end

% Assemble Nodal Masses if Mass is not updated every cycle
if ~isMassUpdateEveryCycle
    massN=zeros(ndofpD,1);
    MMQuad=zeros(ndofpD);
    
    for N=1:NZ
        if ~isRefinedQuad
            localm=OLDdensityZ(N)*GetMassMatrix(refmesh(:,:,N),Basis,massQuadOrder);
        else
            localm=OLDdensityZ(N)*GetMassMatrixR(refmesh(:,:,N),Basis,massQuadOrder);
        end
        massN(Quadmap(:,N))=massN(Quadmap(:,N))+diag(localm,0);
        MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+localm;
    end
end
    

%% === Generation Visualization Phase ===

% Step 2) Now we loop over time steps and solve the hydro equations at each
% Lagrangian time step

% === Lagrangian Phase ===

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

format long
disp(Basis)
disp(['massQuadOrder = ',num2str(massQuadOrder)])
IEinit=sum(massZ.*OLDenergyZ)
KEinit=0.5*(OLDvelocityQuad(1,:)*MMQuad*OLDvelocityQuad(1,:)'+OLDvelocityQuad(2,:)*MMQuad*OLDvelocityQuad(2,:)')
massNinit=sum(MMQuad(:))

% %% === Time Stepping Loop ===
% for cycle=1:maxcycle
%     if stop==true
%         break;
%     end
%     
%     for i=1:NZ
%         len(i)=1/sqrt(ComputeInvLengthSquared(oldmesh(:,:,i)));
%     end
%     
%     % Compute stable time increment
%     if cycle > 1
%         % Ramp up comparison time so that it eventually stops being chosen
%         dt=dtInit;
%         for i=1:NZ
%             zonevel=OLDvelocityQuad(:,topo(:,i));
%             % Maximum velocity in zone
%             umax=sqrt(max(abs(zonevel(1,:)))^2+max(abs(zonevel(2,:)))^2);
%             dt=min([dt,0.1*len(i)/(umax+1e-6)]);
% %             dt=min([dt,1e-3]);
%         end
%     else
%         dt=dtInit;
%     end
%     
%     % Compute sound speed
%     soundSpeedZ=sqrt(gamma*pressureZ./OLDdensityZ);
%     
% %% === Acceleration Phase ===
%     
%     % Compute the hydro corner forces
%     % These are simply the values of the nodal forces at each node of a
%     % given zone
%     if ~isRefinedQuad
%         for i=1:NZ
%             cornerForce(:,:,i)=ComputeCornerForcesFEM(oldmesh(:,:,i),pressureZ(i),GradBasis,stiffQuadOrder);
%         end
%     else
%         for i=1:NZ
%             cornerForce(:,:,i)=ComputeCornerForcesR(oldmesh(:,:,i),pressureZ(i),GradBasis,stiffQuadOrder);
%         end
%     end
%     
%     if hgfrac > 0
%         % Compute the nodal hour glass forces
%         % These are the artificial forces that are designed to "resist" the
%         % hourglass modes
%         % The magnitude of this force can be scaled with the parameter hgfrac
%         for i=1:NZ
%             zonevel=OLDvelocityQuad(:,topo(:,i));
% 
%             thisHGForce=1/(4*2*dt)*massZ(i)*ComputeHGForces(zonevel);
% 
%             % Set the hourglass corner force
%             cornerForce(:,:,i)=cornerForce(:,:,i)+hgfrac*thisHGForce(:,:);
%         end
%     end
%         
%     
%     if Qfrac == 1
%         %Tensor artificial viscosity forces
%         for i=1:NZ
%             zonevel=OLDvelocityQuad(:,Quadmap(:,i));
%             if ~isRefinedQuad
%                 [GradVx(1),GradVx(2)]=LocalInterpolateVelocityGrad2(zonevel(1,:),oldmesh(:,:,i),GradBasis,interpQuadOrder);
%                 [GradVy(1),GradVy(2)]=LocalInterpolateVelocityGrad2(zonevel(2,:),oldmesh(:,:,i),GradBasis,interpQuadOrder);
%                 DivV=GradVx(1)+GradVy(2);
%             else
%                 [dV]=LocalInterpolateVelocityGradR(zonevel,oldmesh(:,:,i),GradBasis,interpQuadOrder);
%                 DivV=dV(1,1)+dV(2,2);
%             end
%             if DivV>=0
%                 thisTenQForce=zeros(ndofpZ,2);
%             else
%                 coef=OLDdensityZ(i)*len(i)*(-qquad*DivV*len(i)+qlin*soundSpeedZ(i));
%                 if ~isRefinedQuad
%                     stiffmat=GetStiffnessMatrix(oldmesh(:,:,i),GradBasis,stiffQuadOrder);
%                 else
%                     stiffmat=GetStiffnessMatrixR(oldmesh(:,:,i),GradBasis,stiffQuadOrder);
%                 end
%                 thisTenQForce=coef*(stiffmat*zonevel');
%             end
%             cornerForce(:,:,i)=cornerForce(:,:,i)-thisTenQForce';
%         end
%     end
%     
%     % Sum-in hydro corner forces to obtain the hydro nodal forces
%     forceN=zeros(2,ndofpD);
%     for i=1:NZ
%         forceN(:,Quadmap(:,i))=forceN(:,Quadmap(:,i))+cornerForce(:,:,i);
%     end
% 
% 
%     if isFullMassMatrixSolve
%         %% Mass Matrix Assembly
%         MMQuad=zeros(ndofpD);
% 
%         for N=1:NZ
%             if ~isRefinedQuad
%                 localm=OLDdensityZ(N)*GetMassMatrix(oldmesh(:,:,N),Basis,massQuadOrder);
%             else
%                 localm=OLDdensityZ(N)*GetMassMatrixR(oldmesh(:,:,N),Basis,massQuadOrder);
%             end
%             MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+localm;
%         end
%         MMX=MMQuad;
%         MMY=MMQuad;
%         
%         % Impose Boundary Conditions
%         % Impose zero boundary at bottom wall
%         for i=1:sum(nbnodes(1))
%             index=bdof(i);
%             MMY(index,:)=0;
%             MMY(:,index)=0;
%             MMY(index,index)=1;
%             forceN(2,index)=0;
%         end
%         % Impose zero boundary at right wall
%         for i=1+nbnodes(1):sum(nbnodes(1:2))
%             index=bdof(i);
%             MMX(index,:)=0;
%             MMX(:,index)=0;
%             MMX(index,index)=1;
%             MMY(index,:)=0;
%             MMY(:,index)=0;
%             MMY(index,index)=1;
%             forceN(:,index)=[0;0];
%         end
%         % Impose zero boundary at top wall
%         for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
%             index=bdof(i);
%             MMX(index,:)=0;
%             MMX(:,index)=0;
%             MMX(index,index)=1;
%             MMY(index,:)=0;
%             MMY(:,index)=0;
%             MMY(index,index)=1;
%             forceN(:,index)=[0;0];
%         end
%         % Impose symmetric boundary at left wall
%         for i=1+sum(nbnodes(1:3)):sum(nbnodes)
%             index=bdof(i);
%             MMX(index,:)=0;
%             MMX(:,index)=0;
%             MMX(index,index)=1;
%             forceN(1,index)=0;
%         end
% 
%         % Solve Linear System
%         accelerationQuad=[(MMX\forceN(1,:)')'; (MMY\forceN(2,:)')'];
%         
%     elseif isMassUpdateEveryCycle
%         % Assemble Nodal Masses
%         massN=zeros(ndofpD,1);
%         MMQuad=zeros(ndofpD);
% 
%         for N=1:NZ
%             if ~isRefinedQuad
%                 localm=OLDdensityZ(N)*GetMassMatrix(refmesh(:,:,N),Basis,massQuadOrder);
%             else
%                 localm=OLDdensityZ(N)*GetMassMatrixR(refmesh(:,:,N),Basis,massQuadOrder);
%             end
%             massN(Quadmap(:,N))=massN(Quadmap(:,N))+diag(localm,0);
%             MMQuad(Quadmap(:,N),Quadmap(:,N))=MMQuad(Quadmap(:,N),Quadmap(:,N))+localm;
%         end
%         
%         % Solve for acceleration
%         for i=1:ndofpD
%             accelerationQuad(:,i)=forceN(:,i)/massN(i);
%         end
%         
%         % Impose Boundary Conditions
%         % Impose symmetric boundary at bottom wall
%         for i=1:sum(nbnodes(1))
%             index=bdof(i);
%             accelerationQuad(2,index)=0;
%         end
%         % Impose zero boundary at right wall
%         for i=1+nbnodes(1):sum(nbnodes(1:2))
%             index=bdof(i);
%             accelerationQuad(:,index)=[0;0];
%         end
%         % Impose zero boundary at top wall
%         for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
%             index=bdof(i);
%             accelerationQuad(:,index)=[0;0];
%         end
%         % Impose symmetric boundary at left wall
%         for i=1+sum(nbnodes(1:3)):sum(nbnodes)
%             index=bdof(i);
%             accelerationQuad(1,index)=0;
%         end
%     else
%         % Solve for acceleration
%         for i=1:ndofpD
%             accelerationQuad(:,i)=forceN(:,i)/massN(i);
%         end
%         
%         % Impose Boundary Conditions
%         % Impose symmetric boundary at bottom wall
%         for i=1:sum(nbnodes(1))
%             index=bdof(i);
%             accelerationQuad(2,index)=0;
%         end
%         % Impose zero boundary at right wall
%         for i=1+nbnodes(1):sum(nbnodes(1:2))
%             index=bdof(i);
%             accelerationQuad(:,index)=[0;0];
%         end
%         % Impose zero boundary at top wall
%         for i=1+sum(nbnodes(1:2)):sum(nbnodes(1:3))
%             index=bdof(i);
%             accelerationQuad(:,index)=[0;0];
%         end
%         % Impose symmetric boundary at left wall
%         for i=1+sum(nbnodes(1:3)):sum(nbnodes)
%             index=bdof(i);
%             accelerationQuad(1,index)=0;
%         end
%     end
% 
%     % Integrate accelerations to get velocities
%     NEWvelocityQuad=OLDvelocityQuad+dt*accelerationQuad;
% 
%     % Integrate velocities to get displacements and update the grid
%    newnodes=oldnodes+dt*NEWvelocityQuad(:,1:NN);
%    newmesh=ComputeMesh(newnodes,topo);
%    
% %% === Work Phase ===   
%    
%    % Calculate the new zone volumes from the new node coordinates
%    for i=1:NZ
%        NEWvolumeZ(i)=ComputeVolume(newmesh(:,:,i));
%        NEWdensityZ(i)=massZ(i)/NEWvolumeZ(i);
%        
%        % Compute the new internal energy
%        % F dot V internal energy update
%        zonevel=0.5*(NEWvelocityQuad(:,Quadmap(:,i))+OLDvelocityQuad(:,Quadmap(:,i)));
%        NEWenergyZ(i)=OLDenergyZ(i)-dt/massZ(i)*sum(dot(zonevel,cornerForce(:,:,i)));
%    end
%    
%    % ===EOS PHASE ===
%    
%    % Evaluate the EOS for an ideal monatomic gas to get the new pressure
%    
%    for i=1:NZ
%        pressureZ(i)=(gamma-1)*NEWdensityZ(i)*NEWenergyZ(i);
%        if pressureZ(i) < 0
%            fprintf('WARNING -- Negative Pressure detected at node %4.0f \n',i)
%        end
%    end
%    
%    % Tally up the energies
%    internalenergy(cycle)=sum(massZ.*NEWenergyZ);
%    if ~isFullMassMatrixSolve
%        kineticenergy(cycle)=0.5*(NEWvelocityQuad(1,:)*MMQuad*NEWvelocityQuad(1,:)'+NEWvelocityQuad(2,:)*MMQuad*NEWvelocityQuad(2,:)');
%        sumMassN(cycle)=sum(massN(:));
%    else
%        kineticenergy(cycle)=0.5*(NEWvelocityQuad(1,:)*MMX*NEWvelocityQuad(1,:)'+NEWvelocityQuad(2,:)*MMY*NEWvelocityQuad(2,:)');
%        sumMassN(cycle)=sum(MMQuad(:));
%    end
%    
%    t=t+dt;
%    if t>=tstop
%        stop=true;
%    end
%    timedata(cycle)=t;
%    
%    % Reset old state variables to new state variables
%    
%    OLDvelocityQuad=NEWvelocityQuad;
%    OLDvolumeZ=NEWvolumeZ;
%    OLDdensityZ=NEWdensityZ;
%    OLDenergyZ=NEWenergyZ;
%    
%    oldmesh=newmesh;
%    oldnodes=newnodes;
%    
%    % === VISUALIZATION AND POST-PROCESSING PHASE ===
%    
% %% =========== PLOT RESULTS ===========
%    if mod(cycle,plotcycle) == 0 || cycle==1
%        fprintf(1,'Cycle: %4.0f,  Time: %8.6f,  Total Energy: %20.18f, Nodal Mass: %20.18f \n',cycle,t,internalenergy(cycle)+kineticenergy(cycle),sumMassN(cycle));
%        plotvarZ=pressureZ;
%        plotvarN=sqrt(NEWvelocityQuad(1,1:NN).^2+NEWvelocityQuad(2,1:NN).^
%        2);
       
%        lineoutZ=plotvarZ(1:ZSTRIDE);
%        lineoutN=plotvarN(1:NSTRIDE);
%        [X,Y,plotvarS]=ZonalToStructured(NZx,NZy,newmesh,plotvarZ);
       
       sfigure(1);
%        subplot(2,2,1)
       trisurf2(topo',newnodes(1,:),newnodes(2,:),plotvarN');
       axis equal tight
       title('Velocity Magnitude','FontWeight','demi','FontSize',14)
       az=0;
       el=90;
       view(az,el);
       colorbar
       
%        subplot(2,2,2)
%        plot(newnodes(1,1:NSTRIDE),lineoutN);
%        title('Velocity Magnitude','FontWeight','demi','FontSize',14)
%        
%        subplot(2,2,3)
%        ZonalSurf(topo',newnodes(1,:),newnodes(2,:),zeros(NN,1),plotvarZ);
%        axis equal tight
%        title('Pressure','FontWeight','demi','FontSize',14)
%        colorbar
% 
%        subplot(2,2,4)
%        hold on
%        plot(timedata(1:cycle),internalenergy(1:cycle)+kineticenergy(1:cycle)-internalenergy(1),'k')
%        plot(timedata(1:cycle),internalenergy(1:cycle)-internalenergy(1),'r')
%        plot(timedata(1:cycle),kineticenergy(1:cycle),'b')
%        xlabel('Time [sec]')
%        ylabel('Energy')
%        legend('Total Energy - Einit','Internal Energy - Einit','Kinetic Energy','Location','NorthOutside')
%      
%        subplot_title(['Initial Disturbance: t = ',num2str(t)],[.1 .9 .8 .05]);
%        
%        if SaveFigures
%            % Save each plot frame for conversion to animated gif
%            saveas(1,['FigureFiles/',figurefile,'/plot',datestr(clock, 'ddHHMMSS'),'.png']);
%        else
%            drawnow
%        end
%    end
   
% end