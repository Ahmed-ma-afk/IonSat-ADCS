function setupCubeSatVisualization()
%setupVisualization Set CubeSat Simulation Visualization Properties

%  Copyright 2019 The MathWorks, Inc.

if ~evalin('base', 'exist(''visSL3D'')')
    evalin('base', 'visSL3D = Simulink.Variant(''variantVisualization == 1'');')
end

if ~evalin('base', 'exist(''visOff'')')
    evalin('base', 'visOff = Simulink.Variant(''variantVisualization == 0'');')
end

visVariant = -1;
if evalin('base', 'exist(''variantVisualization'')')
    switch evalin('base', 'variantVisualization')
        case 0
            visVariant = 2; % No Visualization
        case 1
            visVariant = 1; % Simulink 3D Animation
    end
end

if builtin('license','test','Virtual_Reality_Toolbox')
    [indx,tf] = listdlg('PromptString','Select a Visualization Tool:',...
        'SelectionMode','single', 'ListSize', [200,60],...
        'ListString',{'Simulink 3D Animation','No Visualization'});
    if tf && indx ~= visVariant
        if indx == 1
            assignin('base', 'variantVisualization', 1); % Simulink 3D Animation
        else
            assignin('base', 'variantVisualization', 0); % No Visualization
        end
        if bdIsLoaded('asbCubeSat')
            model_obj = get_param('asbCubeSat','Object');
            model_obj.refreshModelBlocks;
        end
    end
else
    assignin('base', 'variantVisualization', 0);
    if bdIsLoaded('asbCubeSat')
        model_obj = get_param('asbCubeSat','Object');
        model_obj.refreshModelBlocks;
    end
    msgbox('No License found for Simulink 3D Animation.');
end
end
