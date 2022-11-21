%% Clear workspace and figures
clear all
close all
clc

%% Loading the Data
%Loads the necessary physical variables describing the satellite and its
%components
load('SatConstants.mat')
load('data.mat')
load('Simulation.mat')
load('C.mat')
load('F.mat')
load('P.mat')
load('T.mat')
load('workspace.mat')
%Not necessary while perturbations not implemented
load('LOAS.mat')

%% Orbit and attitude parameters
%Orbit: Initialisation of Keplerian parameters
orbit.a = 6678;     %semimajor axis [km]
orbit.e = 0.001;     %eccentricity
orbit.i = 51.6;     %inclination [degrees]
orbit.O = 146;      %Right ascension of the right ascending node [degrees] %max 197, min 300.5 %181
orbit.o = 344;      %Argument of the perigee [degrees]                      %max 90, min 0      %90
orbit.nu = 0;       %True anomaly [degrees]

%Attitude: Initialisation of angles and rotational speeds:
%Initial orientation in ZYX (alpha,beta,gamma) Euler angles [degrees]
att.alpha = 40;     
att.beta = -10;
att.gamma = 60;
att.alpha = 00;     
att.beta = -10;
att.gamma = 00;


%Initial angular velocities in each axis (x,y,z) of body frame [degrees/sec]
att.wx0 = -0.01;        
att.wy0 = 0.01;
att.wz0 = 0.01;

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

%Initialisation of date
date.year = 2022;
date.month = 1;
date.day = 1;
date.hours = 0;
date.minutes = 0;
date.seconds = 0;


%% Declare variables
%time step based on gyro sampling frequency: 
%TimeStep = 0.25;        %fixed-step size in solver, Default time step=0.25
TimeStep = 0.5;        %fixed-step size in solver, Default time step=0.25
Torbit=2*pi*sqrt((orbit.a)^3/(3.986004418E5));
N_orbits = 2;           %number of orbits to be simulated
%Time spent performing the simulation in seconds (one orbit is ~5400 s):
t_sim = N_orbits*Torbit;
%t_sim = 10800;
%OrbitSize = 5400;

%% Other blocks configuration
%Determine the order of approximation for IGRF model of the earth's
%magnetic field nmax (must be equal or less than 13)
nmax = 2;

%Determine whether you are using full IGRF model or dipole approximation
%1 for true, 0 for false
Use_IGRF = 1; 
%Use_IGRF = 1; 
%Determine whether you are using full IGRF model or dipole approximation
%for the Kalman Filter
%1 for true, 0 for false
Use_IGRF_KF = 1; 
%Use_IGRF_KF = 1; 


%Decide whether you'd like to use control or not
%A value of 1 (True) entails that you'd like to simulate a realistic model
%for the control system, and 0 (False) means you assume ideal no control conditions
Enable_control = 0;

%Enable whether the sensors are considering eclipses or not
Enable_eclipse = 1;

%Global time step
global dt
%Sampling times:
Tss=1;         %sampling time for MTM, FSS, CSS could be 1 not 0.5
Ts=1;          %sampling time for Gyro

%Define initial values for Kalman filter parameters
x0 = [0.207740030841190 0.530411062089572 -0.675665322016929 0.467957858597196 1e-3 1e-3 1e-3];
P0 = diag([1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);
What0 = [-8.56493764825016e-05 -0.000707789025389356 -0.000815029912348624];
%Identity matrix with same number of dimensions as P
CovId = eye(6);

dt = TimeStep;
sigma_w = 1e-4;
sigma_b = 1e-5;
sigma_s = 1e-3;
sigma_m = 1e-6;
sigma_sc = 1e-1;

Gk = diag([-1 -1 -1 1 1 1]);
Q1 = ((sigma_w^2)*dt + (sigma_b^2)*(dt^3)/3) .* eye(3);
Q2 = -((sigma_b^2)*(dt^2)/2) .* eye(3);
Q3 = Q2;
Q4 = ((sigma_b^2)*dt) .* eye(3);
Qk = [Q1 Q2; Q3 Q4];
GQG = Gk*Qk*(Gk.');

% Calculate the pseudo inverse matrix
Pinv_RW_repartition

%Aerodynamic Torque Constants
%load('LOAS.mat')
load('IonSat_6U.mat');
IonSataero.T = SNAP_aeromodel.T;
IonSataero.pitch = SNAP_aeromodel.pitch;
IonSataero.roll = SNAP_aeromodel.roll;
IonSataero.av_density_vs_alt = SNAP_aeromodel.av_density_vs_alt;
IonSataero.alt_range = SNAP_aeromodel.alt_range;


%% Open the simulink model
IonSatSimulationC

%% Load the simulink model
load_system("IonSatSimulationC.slx")

%% Run the simulink 
simOut = sim("IonSatSimulationC");