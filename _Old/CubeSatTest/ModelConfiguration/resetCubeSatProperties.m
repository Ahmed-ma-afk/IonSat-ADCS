function resetCubeSatProperties()
%resetCubesatProperties Reset default initial orbital state and attitude,
%gains, and mass properties.

% Initial orbital state is taken from ISS trajectory data for
% 2019/004/12:00:00 UTC (Orbit 2981).  Real-time ISS trajectory data is
% available at:
% https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html

%  Copyright 2019 The MathWorks, Inc.

if ~bdIsLoaded('asbCubeSat') || strcmp(get_param('asbCubeSat','SimulationStatus'), 'stopped')
    %% Reset mass properties
    vehicle.mass = 1;                      % CubeSat mass [kg]
    vehicle.inertia = eye(3);              % CubeSat Moments of Inertia [kg*m^2]
    assignin('base', 'vehicle', vehicle);
    
    %% Reset controller properties
    gains.Kp = 0.00001;                    % Controller Proportional Gain
    gains.Ki = 0.0000000001;               % Controller Integral Gain
    gains.Kd = 0.01;                       % Controller Derivative Gain
    assignin('base', 'gains', gains);
    
    %% Reset default initial orbital state and attitude
      
    mode = 'Keplerian Orbital Elements';   % Input Mode
    date0 = [2019 1 4 12 0 0];
    mjd = mjuliandate(date0);
    sim_t0 = juliandate(date0);            % Simulation start date
    epoch = 2451545.0;                     % Epoch [Julian date] (J2000)
    a = 6786233.13;                        % Semi-major Axis [m]
    ecc = 0.0010537;                       % Eccentricity
    incl = 51.7519;                        % Inclination [deg]
    RAAN = 95.2562;                        % Right ascension of the ascending node [deg]
    argp = 93.4872;                        % Argument of perigee [deg]
    nu  = 302.9234;                        % True anomoly [deg]
    euler = [0 0 0];                       % Euler angles [deg];
    pqr = [0 0 -0.05168];                  % Body angular rates [deg/s];
    dAT = 37;                              % Difference between TAI and UTC [s]
    dUT1 = deltaUT1(mjd,'action','none');  % Difference between UTC and UT1 [s]
    pm = polarMotion(mjd,'action','none')*180/pi; % Polar displacement in x and y axis [deg deg]
    dCIP = deltaCIP(mjd,'action','none')*180/pi;  % Adjustment to CIP in x and y axes [deg deg]
    lod = 0;                               % Excess length of day [s]
    blk = 'asbCubeSat/Edit Initial Orbit and Attitude';
    if ~bdIsLoaded('asbCubeSat')
        open_system('asbCubeSat');
    end
    initCond = computeCubeSatState(blk, mode, sim_t0, epoch, a, ecc, incl, RAAN,...
        argp, nu, 0, 0, 0, euler, pqr, dAT, dUT1, pm, dCIP, lod);
    assignin('base', 'initCond', initCond);
else
     msgbox('Configuration cannot be reset while simulation is running.');
end
end