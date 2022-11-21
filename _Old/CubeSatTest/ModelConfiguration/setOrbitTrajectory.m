function setOrbitTrajectory()
%setOrbitTrajectory Set CubeSat Simulation Trajectory
%Characteristics

%  Copyright 2019 MathWorks, Inc.

if ~bdIsLoaded('asbCubeSat')
    open_system('asbCubeSat');
end

open_system('asbCubeSat/Edit Initial Orbit and Attitude','mask');
end