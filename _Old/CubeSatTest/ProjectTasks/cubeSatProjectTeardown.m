function cubeSatProjectTeardown()
%cubeSatProjectTeardown Teardown CubeSat Simulation Model

%  Copyright 2019 MathWorks, Inc.

% Save current model configuration parameters for next session
stored = load('CubeSatStoredValues.mat');

if evalin('base', 'exist(''vehicle'', ''var'')')
    stored.vehicle = evalin('base', 'vehicle');
end
if evalin('base', 'exist(''gains'', ''var'')')
    stored.gains = evalin('base', 'gains');
end

save(fullfile(fileparts(mfilename('fullpath')), 'CubeSatStoredValues.mat'), '-struct', 'stored');
    
    
