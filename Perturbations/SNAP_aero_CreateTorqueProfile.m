%% This script generates torque profiles for SNAP's aerodynamic model.
%
% -----------------------------------------------------------------------------
%   Copyright (c) 2010-2018 Samir A. Rawashdeh
%   Electrical and Computer Engineering
%   University of Michigan - Dearborn
%  
%   All rights reserved. 
%   
%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met:
%   
%       * Redistributions of source code must retain the above copyright 
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright 
%         notice, this list of conditions and the following disclaimer in 
%         the documentation and/or other materials provided with the distribution
%         
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%   POSSIBILITY OF SUCH DAMAGE.
%  
% ----------------------------------------------------------------------------

clear
%addpath('libaero') %already added in MATLAB installation of SNAP


%% Define "home" volume
Resolution = 0.25; % 0.25cm cube per dot
res = 360;  %home volume dimensions (resolution)
home_volume = zeros(res,res,res) * NaN;

%% Define satellite Shape
%%Original 3U
%[origin, home_volume] = draw_cubesat(30, 0, home_volume, Resolution); 
%Modified IonSat
%length=34 cm 
%width in dots = 23cm  (these two things have to be modified inside the
%function drawcubesatv1)
%heigth in dots = 10cm (modified inside)
[origin, home_volume] = draw_cubesatv1(34, 0, home_volume, Resolution); 
%origin  = [res res res]/2;  % overwrite, numerical rounding caused odd results.
%origin is an output of the function draw_cubesat should be [180 180 180]
%This is the CoG:
CoG  = [181 179 182];  % overwrite, numerical rounding caused odd results.

SNAP_aeromodel.PointCloudModel = home_volume;

plot_volume(home_volume)
drawnow

%% rotate volume
%roll = (0:3:180) * pi/180;
%Elevation = (-90:3:90) * pi/180;
Elevation = (-90:3:90) * pi/180;
%pitch = (0:3:360) * pi/180;
Azimuth = (0:3:360) * pi/180;

%added for IonSat
Cd = 2.4;

T = zeros(length(Azimuth),length(Elevation));
for iEle = 1:length(Elevation)
    for iAzi = 1:length(Azimuth)
        
        DCM_Az = angle2dcm( 0, 0, Azimuth(iAzi), 'XYZ');
        DCM_El = angle2dcm( 0, Elevation(iEle), 0, 'XYZ');
        DCM = DCM_El*DCM_Az;   %first rotates azimuth, then elevation (roll pitch yaw)
        rot_volume = rotate_volume(home_volume, DCM, origin);
        % plot_volume(rot_volume)
        
        % Torque in body frame
        %temp = calc_torque_v1(rot_volume, origin, Cd, Resolution); %before 
        temp = calc_torque_v1(rot_volume, CoG, Cd, Resolution); %after
        
        T(iAzi, iEle) =  sqrt(temp(1)^2+temp(2)^2+temp(3)^2);
        DCM_El1 = angle2dcm( 0, -Elevation(iEle), 0, 'XYZ');
        DCM_Az1 = angle2dcm( 0, 0, -Azimuth(iAzi), 'XYZ');
        DCM1 = DCM_Az1*DCM_El1;   %first rotates -elevation , then -azimuth (roll pitch yaw)
        temp = DCM1 * temp';
        T_b_x (iAzi, iEle)= temp(1);
        T_b_y (iAzi, iEle)= temp(2);
        T_b_z (iAzi, iEle)= temp(3);
        disp(['Elevation: ' num2str(Elevation(iEle)*180/pi) ', Azimuth: '  num2str(Azimuth(iAzi)*180/pi) '/360' ', Torque =' num2str(temp')])
        %           pause
    end
    
    %         plot(T)
    
end
%T(T<1e-8) = 0;

SNAP_aeromodel.T = [T']; %this transpose is because of the simulink
SNAP_aeromodel.T_x = T_b_x';
SNAP_aeromodel.T_y = T_b_y';
SNAP_aeromodel.T_z = T_b_z';

SNAP_aeromodel.Az = Azimuth;
SNAP_aeromodel.El = Elevation;                            
SNAP_aeromodel.alt_range = (200:50:700); % km, altitudes
SNAP_aeromodel.lo_density_vs_alt = [1.78e-10 3.35e-11 8.19e-12 2.34e-12 7.32e-13 2.47e-13 8.98e-14 3.63e-14 1.68e-14 9.14e-15 5.74e-15]; %% averages, Kg/m3
SNAP_aeromodel.av_density_vs_alt = [2.53e-10 6.24e-11 1.95e-11 6.98e-12 2.72e-12 1.13e-12 4.89e-13 2.21e-13 1.04e-13 5.15e-14 2.71e-14]; %% averages, Kg/m3
SNAP_aeromodel.up_density_vs_alt = [3.52e-10 1.06e-10 3.96e-11 1.66e-11 7.55e-12 3.61e-12 1.8e-12 9.25e-13 4.89e-13 2.64e-13 1.47e-13 ]; %% averages, Kg/m3

save('IonSat_6U','SNAP_aeromodel')
%load('IonSat_6U','SNAP_aeromodel') %to plot
%end of script

%the following is just to plot, to have an idea of the value of the torque
vel=7725.84; %circular orbital velocity at 300km
dens=8.19e-12;  %atmosphere density at 300km low
dens=1.95e-11;  %atmosphere density at 300km av

%mesh(SNAP_aeromodel.roll*180/pi, SNAP_aeromodel.pitch*180/pi, SNAP_aeromodel.T)
figure()
mesh(SNAP_aeromodel.Az*180/pi, SNAP_aeromodel.El*180/pi, SNAP_aeromodel.T*(dens*vel^2))
set(gcf,'color','w');
            %set(gca,'XTick',[0:30:90]);
            %set(gca,'YTick',[0:30:180]);
            title('Aerodynamic torque profile')
            ylabel('Elevation Angle (degrees)')
            xlabel('Azimuth Angle (degrees)')
            %zlabel('Pitch Torque Factor (N.m / [Velocity^2 * Air Density)')
            zlabel('Torque Factor (N.m)')