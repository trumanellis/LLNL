%% Stable Time Step
% Compute the length scale
for i=1:NZ
    len(i)=1/sqrt(ComputeInvLengthSquared(oldmesh(:,:,i)));
end

% Compute stable time increment
if cycle > 1
    % Ramp up comparison time so that it eventually stops being chosen
    dt=dtInit*cycle^1.2;
    for i=1:NZ
        zonevel=OLDvelocityQuad(:,topo(:,i));
        % Maximum velocity in zone
        umax=sqrt(max(abs(zonevel(1,:)))^2+max(abs(zonevel(2,:)))^2);
        dt=min([dt,0.3*len(i)/(umax+1e-6),dtMax]);
    end
    dt=min([dt,tstop-t]);
else
    dt=dtInit;
end

% Compute sound speed
soundSpeedZ=sqrt(gamma*pressureZ./OLDdensityZ);