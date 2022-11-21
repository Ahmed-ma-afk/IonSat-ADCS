clear
close

alt = 300;              %Altitude.....................(Km)							
ecc = 0.001;            %Eccentricity											    
inc =  51.6;            %Inclination..................(deg)		
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
