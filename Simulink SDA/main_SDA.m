% MAIN CODE TO BE EXECUTED
clear 

%Orbit simulation
%The object is assumed submitted only to gravity and it's mass is m in kg
m = 12; %mass (kg)
%Initialisation of Keplerian parameters
a = 6678;     %semimajor axis [km]
e = 0.001;    %eccentricity
i = 98;     %inclination [degrees]
O =  10;      %Right ascension of the right ascending node [degrees] %max 197, min 300.5 %181
o = 90;      %Argument of the perigee [degrees]                      %max 90, min 0      %90
nu = 0;       %True anomaly [degrees]
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

%Initialisation of date
date.year = 2024;
date.month = 3;
date.day = 21;
date.hours = 0;
date.minutes = 0;
date.seconds = 0;

%Simulation time parameters
N_orbits = 3;
%Torbit=90*60;       %1 orbit approx. 90 minutes
Torbit=2*pi*sqrt(a^3/(3.986004418E5));
%tsimulation=60*45;   %in [s] 1 orbit
tsimulation=N_orbits*Torbit;
%tsimulation=10000;
%tsimulation=2700;
delta_t = 0.5; %simulation time step (seconds)

%%% MAGNETIC FIELD MODELS - GET PROPER IGRF COEFFICIENTS %%%
date_IGRF = [date.year,date.month,date.day];
%"REAL" Magnetic Field
%Determine whether you are using full IGRF model or dipole approximation
%1 for full model, 0 for dipole approximation
Use_IGRF = 0; 

nmax = 2;

%"EMBEDDED" Magnetic Field
%Determine whether you are using full IGRF model or dipole approximation
%1 for true, 0 for false
Use_IGRF_KF = 0; 
%Determine the order of approximation for embedded IGRF model
nmax_KF = 2;


%Extended Kalman Filter
global dt Tss Ts q_hat0 b_hat0 W6x6 Vmtm Vfss Vcss  
%global PSDgyro PSDmtm PSDfss PSDcss PSDbias
%Sensor parameters
%Sampling times:
Tss=0.5;         %sampling time for MTM, FSS, CSS could be 1 not 0.5
Ts=0.5;          %sampling time for Gyro
%Sensor Noises
%a) Magnetometer (MTM)
sat.sensors.mag_sigma=5.0e-08;
PSDmtm=sat.sensors.mag_sigma^2*Tss;
%b) Gyrometer 
sat.sensors.gyro_sigma=2.620e-04; %original value
sat.sensors.gyro_sigma=0.2e-04; %to test
sat.sensors.gyro_bias=1.160e-04; %original value
sat.sensors.gyro_bias=0.160e-05; %to test
PSDgyro=sat.sensors.gyro_sigma^2*Ts;
PSDbias=sat.sensors.gyro_bias^2*Ts;
b_offset=2e-4;
b_offset=2e-2;
b_offset=sat.sensors.gyro_bias;
%b_offset=0e-2;
b_offset2=-1e-2;    %for the triangular signal
b_offset2=0;    %for the triangular signal
%b_offset2=sat.sensors.gyro_bias;
%c) Fin Sun Sensor (FSS)
sat.sensors.sun_sigma=0.001164;
PSDfss=sat.sensors.sun_sigma^2*Tss;
%d) Coarse Sun Sensor (CSS)
sat.sensors.sun_coarse_sigma=0.1745;
PSDcss=sat.sensors.sun_coarse_sigma^2*Tss;


%Other Kalman Filter Parameters:
%initial values
dt=Ts;
q_hat0=[1;0;0;0]; 
%q_hat0=[0.6;0.8;0.13;0.05]; q_hat0=q_hat0/norm(q_hat0); 
q_hat0=[0.56 0.81 0.10 0.09]; q_hat0=q_hat0/norm(q_hat0); %ok
% q_hat0=[0.56 0.81 0.05 0.09]; q_hat0=q_hat0/norm(q_hat0); %ok ok
% q_hat0=[0.92 0.32 0.05 0.14]; q_hat0=q_hat0/norm(q_hat0); %very close 0
b_hat0=[0;0;0];
%b_hat0=[0.2e-3;0.2e-3;0.2e-3];


%System state noise matrix, the covariance of the process noise matrix Wd
Vg=sat.sensors.gyro_sigma^2;  %ok
Vg=(5e-5)^2;
Vb=sat.sensors.gyro_bias^2;   %ok
Vb=(5.16e-6)^2; %OK
Vb=(0.16e-5)^2; %OK
%Vb=1e-12;

Q11 = ((Vg*dt+(1/3)*Vb*dt^3))*eye(3);
Q12 = (-1/2*Vb*dt^2)*eye(3);
Q22 = (Vb*dt)*eye(3);
Q = [ Q11 , Q12 ; Q12 , Q22];

G=diag([-1 -1 -1 1 1 1]);
W6x6 = G*Q*G';

%for tests
%A) FSS
% Vfss=sat.sensors.sun_sigma^2;
% Vfss=0.005^2; %repeat again here modify %ok 0.0012
% Vfss=0.04^2; %repeat again here modify %ok 0.0012
% Vfss=0.002^2; %repeat again here modify %ok 0.0012
% V=Vfss*eye(3);

%B) MTM and FSS
% Vmtm=sat.sensors.mag_sigma^2; 
% Vmtm=(2e-7)^2;
% Vfss=sat.sensors.sun_sigma^2;
% Vfss=0.002^2; %repeat again here modify %ok 0.0012
% V=[Vmtm*eye(3),zeros(3);zeros(3),Vfss*eye(3)]; 

%C) MTM, CSS and FSS
Vmtm=sat.sensors.mag_sigma^2; 
%Vmtm=(2e-7)^2;
Vfss=sat.sensors.sun_sigma^2;
%Vfss=0.002^2; %repeat again here modify %ok 0.0012
Vcss=sat.sensors.sun_coarse_sigma^2;
Vcss=0.01745^2;
%Vcss=(1e-1)^2;

%% Simulation
MEKFsim_v2R2019a         %MATLAB R2019a
%MEKFsim_v2              %MATLAB R2020b
data=sim('MEKFsim_v2R2019a');
%data=sim('MEKFsim_v2');

%% Run post-processing functions
timesec = data.tout;
sat_Q = getdatasamples(data.sat_Q,(1:length(timesec)));
sat_omega = getdatasamples(data.sat_omega,(1:length(timesec)));
% sat_pos = getdatasamples(simOut.sat_pos,(1:t_sim));
% sat_speed = getdatasamples(simOut.sat_speed,(1:t_sim));
Q_est = getdatasamples(data.Q_est,(1:length(timesec)));
W_est = getdatasamples(data.W_est,(1:length(timesec)));
Q_error = getdatasamples(data.Quat_error,(1:length(timesec)));
Gyro_bias = getdatasamples(data.gyro_bias,(1:length(timesec)));
G_bias=zeros(length(timesec),3);
for i=1:length(timesec)
    G_bias(i,:)=Gyro_bias(:,1,i);
end
b_est = getdatasamples(data.b_est,(1:length(timesec)));
P =  getdatasamples(data.Pcor,(1:length(timesec)));
SunFine = getdatasamples(data.FSS,(1:length(timesec)));
SunCoarse = getdatasamples(data.CSS,(1:length(timesec)));
b_est = getdatasamples(data.b_est,(1:length(timesec)));
tplot=timesec;  %plot in seconds, min, etc...

figure(1)
set(gcf,'color','w');
    subplot(4,1,1)
        plot(tplot,sat_Q(:,1),'b','LineWidth',1)
        hold on;
        plot(tplot,Q_est(:,1),'r','LineWidth',1)
        plot(tplot,Q_est(:,1)+3*(P(:,1)),'m--','LineWidth',1)
        plot(tplot,Q_est(:,1)-3*(P(:,1)),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Q0: Comparison of the estimated and real values for the satellite quaternion')
        ylabel('')
        xlabel('time in seconds')
        grid on
    subplot(4,1,2)
        plot(tplot,sat_Q(:,2),'b','LineWidth',1)
        hold on;
        plot(tplot,Q_est(:,2),'r','LineWidth',1)
        plot(tplot,Q_est(:,2)+3*(P(:,1)),'m--','LineWidth',1)
        plot(tplot,Q_est(:,2)-3*(P(:,1)),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Q1: Comparison of the estimated and real values for the satellite quaternion')
        ylabel('')
        xlabel('time in seconds')
        grid on
    subplot(4,1,3)
        plot(tplot,sat_Q(:,3),'b','LineWidth',1)
        hold on;
        plot(tplot,Q_est(:,3),'r','LineWidth',1)
        plot(tplot,Q_est(:,3)+3*(P(:,2)),'m--','LineWidth',1)
        plot(tplot,Q_est(:,3)-3*(P(:,2)),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Q2: Comparison of the estimated and real values for the satellite quaternion')
        ylabel('')
        xlabel('time in seconds')
        grid on
    subplot(4,1,4)
        plot(tplot,sat_Q(:,4),'b','LineWidth',1)
        hold on;
        plot(tplot,Q_est(:,4),'r','LineWidth',1)
        plot(tplot,Q_est(:,4)+3*(P(:,3)),'m--','LineWidth',1)
        plot(tplot,Q_est(:,4)-3*(P(:,3)),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Q3: Comparison of the estimated and real values for the satellite quaternion')
        ylabel('')
        xlabel('time in seconds')
        grid on
            
figure(3)
set(gcf,'color','w');
    subplot(3,1,1)
        plot(G_bias(:,1),'b','LineWidth',1)
        hold on;
        plot(b_est(:,1),'r','LineWidth',1)
%         plot(b_est(:,1)+3*P(:,4),'m--','LineWidth',1)
%         plot(b_est(:,1)-3*P(:,4),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Comparison of the estimated and real values for the bias in gyro: Bias x')
        ylabel('rad^2/s')
        xlabel('time in seconds')
        grid on
    subplot(3,1,2)
        plot(G_bias(:,2),'b','LineWidth',1)
        hold on;
        plot(b_est(:,2),'r','LineWidth',1)
%         plot(b_est(:,2)+3*P(:,5),'m--','LineWidth',1)
%         plot(b_est(:,2)-3*P(:,5),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Bias y')
        ylabel('rad^2/s')
        xlabel('time in seconds')
        grid on
    subplot(3,1,3)
        plot(G_bias(:,3),'b','LineWidth',1)
        hold on;
        plot(b_est(:,3),'r','LineWidth',1)
%         plot(b_est(:,3)+3*P(:,6),'m--','LineWidth',1)
%         plot(b_est(:,3)-3*P(:,6),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Bias z')
        ylabel('rad^2/s')
        xlabel('time in seconds')
        grid on

figure(4)
set(gcf,'color','w');
    subplot(3,1,1)
        plot(sat_omega(:,1),'b','LineWidth',1)
        hold on;
        plot(W_est(:,1),'r','LineWidth',1)
        plot(W_est(:,1)+3*P(:,4),'m--','LineWidth',1)
        plot(W_est(:,1)-3*P(:,4),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Comparison of the estimated and real values for the satellite angular velocity: Wx')
        ylabel('')
        xlabel('time in seconds')
        grid on
    subplot(3,1,2)
        plot(sat_omega(:,2),'b','LineWidth',1)
        hold on;
        plot(W_est(:,2),'r','LineWidth',1)
        plot(W_est(:,2)+3*P(:,5),'m--','LineWidth',1)
        plot(W_est(:,2)-3*P(:,5),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Wy')
        ylabel('')
        xlabel('time in seconds')
        grid on
    subplot(3,1,3)
        plot(sat_omega(:,3),'b','LineWidth',1)
        hold on;
        plot(W_est(:,3),'r','LineWidth',1)
        plot(W_est(:,3)+3*P(:,6),'m--','LineWidth',1)
        plot(W_est(:,3)-3*P(:,6),'m--','LineWidth',1)
        legend('Real value','Estimated value')
        title('Wz')
        ylabel('')
        xlabel('time in seconds')
        grid on
        
        
figure(5)
set(gcf,'color','w');
    subplot(3,1,1)
        plot(W_est(:,1)-sat_omega(:,1),'b','LineWidth',1)
        hold on;
        plot(3*P(:,4),'m--','LineWidth',1)
        plot(-3*P(:,4),'m--','LineWidth',1)
        legend('Error estimated value')
        title('Error of the angular velocity')
        ylabel('W_x')
        xlabel('time in seconds')
        grid on
    subplot(3,1,2)
        plot(W_est(:,2)-sat_omega(:,2),'b','LineWidth',1)
        hold on;
        plot(3*P(:,5),'m--','LineWidth',1)
        plot(-3*P(:,5),'m--','LineWidth',1)
        legend('Error estimated value')
        title('Error of the angular velocity')
        ylabel('W_y')
        xlabel('time in seconds')
        grid on
        grid on
    subplot(3,1,3)
        plot(W_est(:,3)-sat_omega(:,3),'b','LineWidth',1)
        hold on;
        plot(3*P(:,6),'m--','LineWidth',1)
        plot(-3*P(:,6),'m--','LineWidth',1)
        legend('Error estimated value')
        title('Error of the angular velocity')
        ylabel('W_z')
        xlabel('time in seconds')
        grid on
        grid on
        
 figure(6)
 set(gcf,'color','w');
    subplot(2,1,1)
    %plot(tplot, rad2deg(quat2eul(Q_error)),'LineWidth',1)
    error=rad2deg(2*acos(Q_error(:,1)));
    for i=1:length(error)
        if error(i)>180
            error(i)=360-error(i);
        end
    end
    plot(tplot, rad2deg(2*acos(Q_error(:,1))),'LineWidth',1)
    title('Error angles (current implementation)')
    ylabel('Degrees')
    xlabel('Time')
    ylim([-10,10])
    yline(5, 'r--', 'LineWidth', 1);
    yline(-5, 'r--', 'LineWidth', 1);
    legend('Error','Requirement error +','Requirement error -')
    grid on
    subplot(2,1,2)
    plot(tplot, SunCoarse,'m--','LineWidth',2)
    hold on
    plot(tplot, SunFine,'LineWidth',2)
    title('Error angles (current implementation)')
    ylabel('Degrees')
    xlabel('Time')
    ylim([-0.5,1.5])
    legend('Coarse sun sensor active','Fine sun sensor active')
    grid on