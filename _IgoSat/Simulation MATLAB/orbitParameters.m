% ------------------------------------------------
%   PROJET  - Attitude control of IGOSAT 
%   Antoine MARTIN-LAGARDE / Oct 2018 - Jan 2019       
% ------------------------------------------------
%   Author : Antoine MARTIN-LAGARDE
%   Date   : 14/12/2018
% ------------------------------------------------

clear; close all; clc;


load IGRFsimCoefs
load IGRFembCoefs


addpath MySimulinkLib
slblocks;


% Constants
%------------------------------------------------------------
mu = 3.98600436e5;  % Earth attraction constant (km^3/sec^-2)
eqRad = 6371*1e3;  % Mean earth radius (m)
deg2rad = pi/180;
rad2deg = 180/pi;

% Initialise Keplar elements for the orbit
%-------------------------------------------------------------
semimajoraxis = 650.0 + eqRad / 1000;   % km
eccentricity  = 0;
inc_SunSynchronous = 97.78739;  % deg
inclinationOrbit   = 90;        % deg     98.602802
w_init = 0;                     % argument of perigee(rad)
Omega_init = 0;                 % longitude of the ascending node(rad)
M = 0;                          % mean anamoly (rad)
t_init= 0;                      % time tle

kep_init = [ semimajoraxis;        % (km)
    inclinationOrbit*deg2rad;              % rad
    Omega_init;                            % rad
    w_init;                                % rad
    eccentricity;
    M;                                     % rad
    t_init];                               % s

omega_o = sqrt(mu/abs(semimajoraxis)^3);      % Orbit rate

% Filtrage

tho=5;
b0=2/11;
b1=-2/11;
a1=-9/11; 


% Time parameters
%--------------------------------------------------------------------------
tme_init = [2018 06 07 16 15 0];                  % time of launch in UTC
tstart = Dt2jD(tme_init);                       % start time in Julian days
tstep = 10;                                     % sampling time
torbit = 2*pi/(sqrt(mu/abs(kep_init(1))^3));    % in sec
torbit_hr = torbit/3600;                        % in hrs
tfinal = 30*torbit;                              % time for 30 orbits(sec)
% timeD = tstart + [0:tstep/86400:tfinal/86400];  % Time grid in Julian days
% timeS = 0:tstep:tfinal;                         % Time grid in seconds
% D_jDdatenum=Dt2jD([0 1 2 0 0 0]);

% satellite characteristics
%--------------------------------------------------------------------------

% mass=3.77;                                      % in kg
% lenth=34.5*10^-2;                               % in m
% width=10*10^-2;                                 % in m
% Iz =mass/6*(width^2);
% Ix=mass*(lenth^2+width^2)/12;
% Iy=Ix;
% I=diag([Ix,Iy,Iz]);
% Ixy=Ix/Iy;
I=10^-6*[29290.7 53.8 0.5; 53.8 28639.4 0; 0.5 0 5490.1];
% spinRate=0.1;

% Initial conditions
Roll=0;
Pitch=90;
Yaw=0;
wx=30;
wy=0;
wz=0;
nIGRF=4;

% Reaction wheel
Hw = [0; 1.77*1e-3; 0];   % Momentum storage at max rpm --> Rotation axe on y satellite
Iw = Hw * 180 / (8000*pi);
% Iw=[0; 2.2125*1e-7; 0];         % Reaction wheel inertial matrix
% Actuator characteristics
%--------------------------------------------------------------------------
% Actuator
mr_max = 0.2;               % maximum magnetic moment of torque rod (A-m^2)
mc_max = 0.2;               % maximum magnetic moment of air-coil (A-m^2)
% i_max = 250*1e-3;            % maximum current (Amps)
vr_max=2.5;                  % maximum Voltage to torque rod (V)
vc_max=5;                    % maximum voltage to air-coil (V)

% Area of each coil
% Ax = 60*10*10^-6; % m^2
% Ay = 60*10*10^-6; % m^2
% Az = 90^2*10^-6;  % m^2

% Resistance in each coil
Rx = 30;
Ry = 30;
Rz = 83;

% Gain de control
kd = 2e-2*min(diag(I))*0.2/3e-9;
Kp=3.29838*1e-8*0.1/1e-12;
Kd=3.68819*1e-4*0.1/8e-11;

% K matrix ?? controll allocation matrix for regulator
% K = [(Nx*Ax)/Rx 0 0; 0 (Ny*Ay)/Ry 0; 0 0 (Nz*Az)/Rz];
% inK=inv(K);

% residual magnetic moment
%--------------------------------------------------------------------------
mr = [0.001; -0.001; 0.005];
% bavg = 1.0e-05 *[-0.0067;0.3398;0.0344];

% Sampling time
%--------------------------------------------------------------------------
% st=1;   % sensor noise 
% Ts=1;   % sampling time of controller & Attitude estimation
ss=0.1;   % Simulation step
cr = 1;   % Control rate

% Sensor & Actuator noise
%--------------------------------------------------------------------------
mn = 3e-18;   % Magnetometer noise (Tesla) : PSD
% sn =1e-8;    % sun sensor noise
% en =1e-8 ;   % earth sensor noise
% gn = 1e-4*[0.001 0.001 0.001]; % gyro noise(rad/sec)
% mtn=1e-9;    % magnetorquer noise