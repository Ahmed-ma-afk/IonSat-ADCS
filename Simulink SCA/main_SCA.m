%% Initialization parameters
clear
close						

%Orbit: Initialisation of Keplerian parameters
orbit.a = 6678;     %semimajor axis [km]
orbit.e = 0.001;     %eccentricity
orbit.i = 51.6;     %inclination [degrees]
orbit.O = 146;      %Right ascension of the right ascending node [degrees] %max 197, min 300.5 %181
%orbit.O = 91;      %Right ascension of the right ascending node [degrees] %max 197, min 300.5 %181
orbit.o = 344;      %Argument of the perigee [degrees]                      %max 90, min 0      %90
%orbit.o = 90;      %Argument of the perigee [degrees]                      %max 90, min 0      %90
orbit.nu = 0;       %True anomaly [degrees]


%Attitude: Initialisation of angles and rotational speeds:
%Initial orientation in ZYX (alpha,beta,gamma) Euler angles [degrees]
att.alpha = 40;     
att.beta = -10;
att.gamma = 60;

%to test large error
att.alpha = 140;     
att.beta = -70;
att.gamma = 90;
%Initial angular velocities in each axis (x,y,z) of body frame [degrees/sec]
att.wx0 = -4;        
att.wy0 = 5;
att.wz0 = 3;

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

%Load other parameters
load('SatConstants.mat')
load('workspace.mat')

%% Other blocks configuration
%Determine the order of approximation for IGRF model of the earth's
%magnetic field nmax (must be equal or less than 13)
nmax = 2;

%Determine whether you are using full IGRF model or dipole approximation
%1 for true, 0 for false
Use_IGRF = 0; 

% Calculate the pseudo inverse matrix
Pinv_RW_repartition

%ControlR2019A
Control_v2

