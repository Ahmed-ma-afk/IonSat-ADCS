clear
close

alt = 300;              %Altitude.....................(Km)							
ecc = 0.001;            %Eccentricity											    
inc =  98;              %Inclination..................(deg)		
RAAN = 30;              %Right Asc. of Ascending Node.(deg)	
w = 25;                 %Argument of perigee..........(deg)	
nu = 75;                %Satellite position...........(deg)							

inc = rad2deg(inc);
RAAN = rad2deg(RAAN);
w = rad2deg(w);
nu = rad2deg(nu);

%time
TimeStep = 1;        %fixed-step size in solver, Default time step=0.25
%TimeStep = 0.5;        %fixed-step size in solver, Default time step=0.25
Torbit=2*pi*sqrt((alt+6378.1)^3/(3.986004418E5));
N_orbits = 2;           %number of orbits to be simulated
%Time spent performing the simulation in seconds (one orbit is ~5400 s):
t_sim = N_orbits*Torbit;

%Initialisation of date
date.year = 2022;
date.month = 3;
date.day = 21;
date.hours = 0;
date.minutes = 0;
date.seconds = 0;

%Load Satellite Constants
load('SatConstants.mat')

%%% MAGNETIC FIELD MODELS - GET PROPER IGRF COEFFICIENTS %%%
date_IGRF = [date.year,date.month,date.day];
%"REAL" Magnetic Field
%Determine whether you are using full IGRF model or dipole approximation
%1 for full model, 0 for dipole approximation
Use_IGRF = 0; 
nmax = 2;

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
%Solar radiation pressure torque
F_s = 1367;     %solar constant in [W/m^2]
c = 299792458;  %light speed in [m/s]
q = 0.6;        %reflectance factor 
%q = (400/((alt+6378)/6378)^2)/F_s;        %reflectance factor 
Cd = 2.4;


distmodel_v5