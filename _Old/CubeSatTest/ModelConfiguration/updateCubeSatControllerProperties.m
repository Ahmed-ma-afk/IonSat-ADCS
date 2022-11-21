function updateCubeSatControllerProperties()
% updateCubeSatControllerProperties Set CubeSat Simulation Model Controller
% Properties

%  Copyright 2019 MathWorks, Inc.

if evalin('base', 'exist(''gains'')') && evalin('base', 'isfield(gains, ''Kp'')') ...
        && evalin('base', 'isfield(gains, ''Ki'')') && evalin('base', 'isfield(gains, ''Kd'')')
    initKp = num2str(evalin('base', 'gains.Kp'));
    initKi = num2str(evalin('base', 'gains.Ki'));
    initKd = num2str(evalin('base', 'gains.Kd'));
else
    stored = load('CubeSatStoredValues.mat', 'gains');
    initKp = num2str(stored.gains.Kp);
    initKi = num2str(stored.gains.Ki);
    initKd = num2str(stored.gains.Kd);
end

prompt = {'CubeSat Proportional Controller Gain (Kp):','CubeSat Integral Controller Gain (Ki):',...
    'CubeSat Derivative Controller Gain (Kd):'};
title = 'CubeSat Controller Properties';
dims = [1 50];
definput = {initKp,initKi,initKd};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
userInput = inputdlg(prompt,title,dims,definput, opts);

if ~isempty(userInput)
    gains.Kp = eval(userInput{1});
    gains.Ki = eval(userInput{2});
    gains.Kd = eval(userInput{3});
    assignin('base', 'gains', gains);
end

end