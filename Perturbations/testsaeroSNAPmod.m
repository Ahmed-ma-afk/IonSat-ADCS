%% create volume & Define "home" volume
clear
Resolution = 0.25; % 0.25cm cube per dot
res = 360;  %home volume dimensions (resolution)
home_volume = zeros(res,res,res) * NaN;

%Define satellite Shape
%length=34 cm 
%width in dots = 23cm  (these two things have to be modified inside the %function drawcubesatv1)
%heigth in dots = 10cm (modified inside)
[origin, home_volume] = draw_cubesatv1(34, 0, home_volume, Resolution); 
%origin is an output of the function draw_cubesat should be [180 180 180]
%This is the CoG:
CoG  = [181 179 182];  % overwrite, numerical rounding caused odd results.
SNAP_aeromodel.PointCloudModel = home_volume;
figure(1)
plot_volume(home_volume)
drawnow

%% Rotation of volume 
Az=150*pi/180;
El=-90*pi/180;

        DCM_Az = angle2dcm( 0, 0, Az, 'XYZ');
        DCM_El = angle2dcm( 0, El, 0, 'XYZ');
        DCM = DCM_El*DCM_Az;   %first rotates azimuth, then elevation (roll pitch yaw)
        rot_volume = rotate_volume(home_volume, DCM, origin);
        
figure(2)
plot_volume(rot_volume)
drawnow

%% calculate torque for some orientations
%Elevation = (-90:3:90) * pi/180;
Elevation = 90 * pi/180;
Azimuth = (0:5:360) * pi/180;

Cd = 2.4;

T = zeros(length(Azimuth),length(Elevation));
for iEle = 1:length(Elevation)
    for iAzi = 1:length(Azimuth)
        
        DCM_Az = angle2dcm( 0, 0, Azimuth(iAzi), 'XYZ');
        DCM_El = angle2dcm( 0, Elevation(iEle), 0, 'XYZ');
        DCM = DCM_El*DCM_Az;   %first rotates azimuth, then elevation (roll pitch yaw)
        rot_volume = rotate_volume(home_volume, DCM, origin);
        figure(2)
        plot_volume(rot_volume)
        
        % Torque in body frame
        temp = calc_torque_v1(rot_volume, CoG, Cd, Resolution);
        
        T(iAzi, iEle) =  sqrt(temp(1)^2+temp(2)^2+temp(3)^2);
        DCM_El1 = angle2dcm( 0, -Elevation(iEle), 0, 'XYZ');
        DCM_Az1 = angle2dcm( 0, 0, -Azimuth(iAzi), 'XYZ');
        DCM1 = DCM_Az1*DCM_El1;   %first rotates -elevation , then -azimuth (roll pitch yaw)
        temp = DCM1 * temp';
        disp(['Elevation: ' num2str(Elevation(iEle)*180/pi) ', Azimuth: '  num2str(Azimuth(iAzi)*180/pi) '/360' ', Torque =' num2str(temp')])
        %           pause
    end
end

normT = [T']; %this transpose is because of the simulink
Azi = Azimuth;
Ele = Elevation;                            


vel=7725.84; %circular orbital velocity at 300km
dens=8.19e-12;  %atmosphere density at 300km

figure(3)
plot(Azi*180/pi,normT*(dens*vel^2))
%mesh(SNAP_aeromodel.Az*180/pi, SNAP_aeromodel.El*180/pi, SNAP_aeromodel.T*(dens*vel^2))
            grid on            
            title('Aerodynamic torque profile')
            %ylabel('Elevation Angle (degrees)')
            xlabel('Azimuth Angle (degrees)')
            ylabel('Torque (N.m)')