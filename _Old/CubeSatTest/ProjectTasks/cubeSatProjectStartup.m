% cubeSatProjectStartup Initialize CubeSat Simulation Model
% Initial orbital state is taken from ISS trajectory data for
% 2019/004/12:00:00 UTC (Orbit 2981).  Real-time ISS trajectory data is
% available at:
% https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html

%  Copyright 2019 The MathWorks, Inc.

CubeSatTimeStep = 1;

%% Load model bus definitions
cubeSatBusDefinitions;

%% Visualization control
visSL3D = Simulink.Variant('variantVisualization == 1');
visOff = Simulink.Variant('variantVisualization == 0');

if builtin('license','test','Virtual_Reality_Toolbox')
    variantVisualization = 1; % Simulink 3D Animation
else
    variantVisualization = 0; % No Visualization
end

%% Load mass properties and controller gains 
load('CubeSatStoredValues.mat');
