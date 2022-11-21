function openVisualizationSettings()
%openVisualizationSettings Change CubeSat Animation Settings

%  Copyright 2019 MathWorks, Inc.

if ~bdIsLoaded('asbCubeSat')
    open_system('asbCubeSat');
end

open_system('asbCubeSat/Visualization','mask');
end