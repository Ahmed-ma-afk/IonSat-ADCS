% cubeSatBusDefinitions Define bus structures used in the CubeSat
% Simulation

%  Copyright 2019 The MathWorks, Inc.

%% Bus object: StatesOutBus

busElemNames = {'V_ecef', 'X_ecef', 'q_ecef2b', 'q_eci2b', 'Euler', 'LatLonAlt',...
    'x_sun_eci'};
busElemDims = {3, 3, 4, 4, 3, 3, 3};

for idx = numel(busElemNames):-1:1
    busElems(idx) = Simulink.BusElement;
    busElems(idx).Name = busElemNames{idx};
    busElems(idx).Dimensions = busElemDims{idx};
    busElems(idx).DimensionsMode = 'Fixed';
    busElems(idx).DataType = 'double';
    busElems(idx).SampleTime = -1;
    busElems(idx).Complexity = 'real';
    busElems(idx).SamplingMode = 'Sample based';
    busElems(idx).Min = [];
    busElems(idx).Max = [];
    busElems(idx).DocUnits = '';
    busElems(idx).Description = '';
end

StatesOutBus = Simulink.Bus;
StatesOutBus.HeaderFile = '';
StatesOutBus.Description = '';
StatesOutBus.DataScope = 'Auto';
StatesOutBus.Alignment = -1;
StatesOutBus.Elements = busElems;

clear busElems busElemNames busElemDims idx;

%% Bus object: ACSOutBus

busElemNames = {'AttitudeError', 'TorqueCmds', 'NavMode'};
busElemType = {'AttitudeErrorBus', 'double', 'double'};
busElemDims = {1, 3, 1};

for idx = numel(busElemNames):-1:1
    busElems(idx) = Simulink.BusElement;
    busElems(idx).Name = busElemNames{idx};
    busElems(idx).Dimensions = busElemDims{idx};
    busElems(idx).DimensionsMode = 'Fixed';
    busElems(idx).DataType = busElemType{idx};
    busElems(idx).SampleTime = -1;
    busElems(idx).Complexity = 'real';
    busElems(idx).SamplingMode = 'Sample based';
    busElems(idx).Min = [];
    busElems(idx).Max = [];
    busElems(idx).DocUnits = '';
    busElems(idx).Description = '';
end

ACSOutBus = Simulink.Bus;
ACSOutBus.HeaderFile = '';
ACSOutBus.Description = '';
ACSOutBus.DataScope = 'Auto';
ACSOutBus.Alignment = -1;
ACSOutBus.Elements = busElems;

clear busElems busElemNames busElemType busElemDims idx;

%% Bus object: FCSOutBus

busElems(1) = Simulink.BusElement;
busElems(1).Name = 'ACSOut';
busElems(1).Dimensions = 1;
busElems(1).DimensionsMode = 'Fixed';
busElems(1).DataType = 'ACSOutBus';
busElems(1).SampleTime = -1;
busElems(1).Complexity = 'real';
busElems(1).SamplingMode = 'Sample based';
busElems(1).Min = [];
busElems(1).Max = [];
busElems(1).DocUnits = '';
busElems(1).Description = '';

FCSOutBus = Simulink.Bus;
FCSOutBus.HeaderFile = '';
FCSOutBus.Description = '';
FCSOutBus.DataScope = 'Auto';
FCSOutBus.Alignment = -1;
FCSOutBus.Elements = busElems;

clear busElems;

%% Bus object: AttitudeErrorBus

busElemNames = {'Roll', 'Pitch', 'Yaw'};

for idx = numel(busElemNames):-1:1
    busElems(idx) = Simulink.BusElement;
    busElems(idx).Name = busElemNames{idx};
    busElems(idx).Dimensions = 1;
    busElems(idx).DimensionsMode = 'Fixed';
    busElems(idx).DataType = 'double';
    busElems(idx).SampleTime = -1;
    busElems(idx).Complexity = 'real';
    busElems(idx).SamplingMode = 'Sample based';
    busElems(idx).Min = [];
    busElems(idx).Max = [];
    busElems(idx).DocUnits = '';
    busElems(idx).Description = '';
end

AttitudeErrorBus = Simulink.Bus;
AttitudeErrorBus.HeaderFile = '';
AttitudeErrorBus.Description = '';
AttitudeErrorBus.DataScope = 'Auto';
AttitudeErrorBus.Alignment = -1;
AttitudeErrorBus.Elements = busElems;

clear busElems busElemNames idx;

%% Bus object: EnvBus

busElemNames = {'envForces_body', 'envTorques_body', 'x_sun_eci', 'earthAngRate'};
busElemDims = {3, 3, 3, 1};

for idx = numel(busElemNames):-1:1
    busElems(idx) = Simulink.BusElement;
    busElems(idx).Name = busElemNames{idx};
    busElems(idx).Dimensions = busElemDims{idx};
    busElems(idx).DimensionsMode = 'Fixed';
    busElems(idx).DataType = 'double';
    busElems(idx).SampleTime = -1;
    busElems(idx).Complexity = 'real';
    busElems(idx).SamplingMode = 'Sample based';
    busElems(idx).Min = [];
    busElems(idx).Max = [];
    busElems(idx).DocUnits = '';
    busElems(idx).Description = '';
end

EnvBus = Simulink.Bus;
EnvBus.HeaderFile = '';
EnvBus.Description = '';
EnvBus.DataScope = 'Auto';
EnvBus.Alignment = -1;
EnvBus.Elements = busElems;

clear busElems busElemNames busElemDims idx;