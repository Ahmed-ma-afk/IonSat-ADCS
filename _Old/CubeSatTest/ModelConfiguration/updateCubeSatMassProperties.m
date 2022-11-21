function updateCubeSatMassProperties()
%updateCubeSatMassProperties Set CubeSat Simulation Model Mass Properties

%  Copyright 2019 MathWorks, Inc.

if evalin('base', 'exist(''vehicle'')') && evalin('base', 'isfield(vehicle, ''mass'')') ...
        && evalin('base', 'isfield(vehicle, ''inertia'')')
    initMass = num2str(evalin('base', 'vehicle.mass'));
    moi = num2str(evalin('base', 'vehicle.inertia'));
else
    stored = load('CubeSatStoredValues.mat', 'vehicle');
    initMass = num2str(stored.vehicle.mass);
    moi = num2str(stored.vehicle.inertia);
end
initMOI = regexprep(['[' moi(1,:) ';' moi(2,:) ';' moi(3,:) ']'],' +',' ');

prompt = {'CubeSat Mass (kg):','CubeSat Moments of Inertia:'};
title = 'CubeSat Mass Properties';
dims = [1 50];
definput = {initMass,initMOI};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
userInput = inputdlg(prompt,title,dims,definput,opts);

if ~isempty(userInput)
    vehicle.mass = eval(userInput{1});
    vehicle.inertia = eval(userInput{2});
    assignin('base', 'vehicle', vehicle);
end
end