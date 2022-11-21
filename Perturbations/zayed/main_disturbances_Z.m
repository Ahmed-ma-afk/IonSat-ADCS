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

%......
date_IGRF = [date.year,date.month,date.day];
Use_IGRF = 0; 
nmax = 2;
%....

%Satelite dipole resudiale
RSD=[0.08,0.08,0.08];
%about the new model
AX=0.022;
AY=0.034;
AZ=0.2108;
OX=0.17;
OY=0.11;
OZ=0.05;
Cdrag=2.6;
%gravity center
OGX=2.5*1e-3;
OGY=-1*1e-2;
OGZ=-2*1e-2;
G=[OGX,OGY,OGZ];
%about deployable_solar_plan
CRS=[0,0.195,0.05];%%the_center_of_the_right_one.
CLS=[0,-0.195,0.05];%%the_center_of_the_left_one.
ESS=0.0578;%%surface_of_one
AZES=0.1156;%surface_of_the_2
%about_the_surface_facing Y
Cmy=[0,-0.11,0];%the_center_of_the_one_facing_-Y.

load('LOAS.mat')
load('SatConstants.mat')

distmodelR2019a_Zayed
%distmodelR2019_Zayed