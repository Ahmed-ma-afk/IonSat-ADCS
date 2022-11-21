%% Initialization parameters
clear
close						

%Orbit: Initialisation of Keplerian parameters
orbit.a = 6678;     %semimajor axis [km]
orbit.e = 0.001;    %eccentricity
orbit.i = 98;     %inclination [degrees]
orbit.O = 10;      %Right ascension of the right ascending node [degrees] %max 197, min 300.5 %181
orbit.o = 90;      %Argument of the perigee [degrees]                      %max 90, min 0      %90
orbit.nu = 0;       %True anomaly [degrees]
%beta angle  ~9 deg. 21/3, i=98, O=10, 
%beta angle ~29 deg. 21/3, i=98, O=30,
%beta angle ~48 deg. 21/3, i=98, O=50,  
%i = 98: notice that for the date 21/3 (equinox) the RAAN (O) is close to the beta angle.

%POINTING MODE
MODE = 2;   %"sun-aero" pointing mode
% 1: "orbital" Reference quaternion is aligned with ORF. 
% 2: "sun-aero" Reference quaternion is such that x is aligned with velocity 
%and z is aligned as best as possible with the sun direction to maximize the power generation 
% 3: "sun pointing" Reference quaternion is such that z is aligned with with 
%the sun direction to maximize the power generation and x is aligned as best 
%as possible with the velocity direction. 
% 4: "aero-drag" Reference quaternion is similar than in Case 1, but rotated 
%90° along the y axis, therefore, the reference quaternion is an attitude for 
%maximizing the drag surface. 
% 5: "retrogade firing" Reference quaternion is similar than in Case 1, but 
%rotated 180° along the z axis, therefore, the reference quaternion is an 
%attitude for retrograde propulsion. (??should be like case 2 but rotated?)
% 6: Reference quaternion is static [1 0 0 0], but this time is the only case 
%where the state of the B-dot is enabled (1). 
% 7: Reference quaternion is a custom quaternion that has to be defined by the 
%user as an input to the mission block. 

%Attitude: Initialisation of angles and rotational speeds:
% %Initial orientation in ZYX (alpha,beta,gamma) Euler angles [degrees]
% att.alpha = 40;     
% att.beta = -10;
% att.gamma = 60;

% %to test large error
% att.alpha = 140;     
% att.beta = -70;
% att.gamma = 90;
% %Initial angular velocities in each axis (x,y,z) of body frame [degrees/sec]
% att.wx0 = -4;        
% att.wy0 = 5;
% att.wz0 = 3;

%to test other errors
att.alpha = 50;         %Initial orientation Yaw [deg]
att.beta = 20;         %Initial orientation Pitch [deg]
att.gamma = 160;         %Initial orientation Roll [deg]
%Initial angular velocities in each axis (x,y,z) of body frame [degrees/sec]
att.wx0 = -0.1;        
att.wy0 = 0.1;
att.wz0 = 0.1;
%works: 10, -20, 30 and -1, 1, 1. mode 2, not 4

%time
TimeStep = 1;        %fixed-step size in solver, Default time step=0.25
Torbit=2*pi*sqrt((orbit.a)^3/(3.986004418E5));
N_orbits = 3;           %number of orbits to be simulated
%Time spent performing the simulation in seconds (one orbit is ~5400 s):
t_sim = N_orbits*Torbit;

%Initialisation of date
%date.year = 2022;
date.year = 2024;
%date.month = 1;
date.month = 3;
%date.day = 1;
date.day = 21;
date.hours = 0;
date.minutes = 0;
date.seconds = 0;

%% Load other parameters
load('SatConstants.mat')
%load('workspace.mat') %for 2021 simulation
sat.inertia = [0.06   0   0;...   %already saved
                0   0.09  0;...
                0     0   0.14];
sat.mass = 12;      %Satellite Mass [kg]
sat.CoG = [0;5.1;4.2]/1000;   %Satellite Center of Gravity [m] in the BRF 
% Other blocks configuration

%%% GET PROPER IGRF COEFFICIENTS %%%
date_IGRF = [date.year,date.month,date.day];
%Determine the order of approximation for IGRF model of the earth's
%magnetic field nmax (must be equal or less than 13)
nmax = 2;
%Determine whether you are using full IGRF model or dipole approximation
%1 for full model, 0 for dipole approximation
Use_IGRF = 0; 


%RWA
% Calculate the pseudo inverse matrix
Pinv_RW_repartition     %executes the following code lines:
% alpha = 26.6 * pi/180;      %from deg to rad
% sat.wheel.repartition_matrix_4RW = [[cos(alpha),      0,    -cos(alpha),      0    ];...
%                                     [     0,     cos(alpha),     0,     -cos(alpha)];...
%                                     [sin(alpha), sin(alpha), sin(alpha), sin(alpha)]];
% sat.wheel.repartition_metrix_4RW_inverse = pinv(sat.wheel.repartition_matrix_4RW);

%% Other disturbance block parameters
%%% Model of Perturbation Torques
%Magnetic Torque Residual dipole
sat.residual_dipole=[0.08;0.08;0.08]; %residual magnetic dipole in [A*m^2]
%Aerodynamic Torque Constants
%load('LOAS.mat')`
load('IonSat_6U.mat');
IonSataero.T = SNAP_aeromodel.T;
IonSataero.El = SNAP_aeromodel.Az;
IonSataero.Az = SNAP_aeromodel.El;
IonSataero.Tx = SNAP_aeromodel.T_x;
IonSataero.Ty = SNAP_aeromodel.T_y;
IonSataero.Tz = SNAP_aeromodel.T_z;
IonSataero.av_density_vs_alt = SNAP_aeromodel.av_density_vs_alt;
IonSataero.alt_range = SNAP_aeromodel.alt_range;

%Thruster direction and activation:
sat.thruster.force=0.00075; % Force of thruster in [N]
d=1*pi/180;           %deviation in rads of the thruster wrt to x axis in plane XY
alpha=80*pi/180;       %rotation of deviation around x axis
    R=[ 1    0            0;...
        0    cos(alpha)  -sin(alpha);...
        0    sin(alpha)   cos(alpha)];
sat.thruster.dir=R*[cos(d);sin(d);0];   %Direction of thrust (unit vector)
%Thruster applied point in BRF:
x_force_thruster=182.95;       %Distance from the geometric center x axis in [mm]
y_force_thruster=0;      %Distance from the geometric center y axis in [mm]
z_force_thruster=0;      %Distance from the geometric center z axis in [mm]
sat.thruster.point=[x_force_thruster;y_force_thruster;z_force_thruster]/1e3;    
                            %Distance from the geometric center (BRF) in [m]
%other parameters for thruster firing
sat.thruster.duration=3000;         %duration of the thrust in [seconds]
sat.thruster.wait=7000;             %waiting time between thrusts in [s]
sat.thruster.firstimpulse=4000;     %first thrust after start sim in [s]
%% Simulation
%Control_v3R2019a           %MATLAB R2019a
%Control_v4                 %MATLAB R2020b
Control_v4_desat0           %MATLAB R2020b

