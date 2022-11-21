clear

d=pi/180;
Time= 50 * 60; %50 minutes

d_alpha=3*pi/180;
a=0:d_alpha:2*pi;

x=0;
y=0.00381;
z=0.00823;

% y=0.0041;
% z=0.0085;
% 
% a=0;

MMSx=(cos(a)*sin(d)*z-sin(a)*sin(d)*y)*Time*0.75;
MMSy=(sin(a)*sin(d)*0.171-cos(d)*z)*Time*0.75;
MMSz=(cos(d)*y-cos(a)*sin(d)*0.171)*Time*0.75;

figure()
set(gcf,'color','w');
    subplot(3,1,1)
    plot(a*180/pi,MMSx);
    grid on
    subplot(3,1,2)
    plot(a*180/pi,MMSy);
    grid on
    subplot(3,1,3)
    plot(a*180/pi,MMSz);
    grid on
    
%% Calc 2
clear

d=pi/180;           %deviation in rads of the thruster wrt to x axis in plane XY
Time= 50 * 60;      %time in seconds (50 minutes)

%the deviation is of d=1 degree
F=0.75*[cos(d);sin(d);0]*1e-3;
y=3.81;             %Distance from the geometric center y axis in mm
z=8.23;             %Distance from the geometric center z axis in mm
D=[171;y;z]/1e3;    %Distance from the geometric center in m

%the rotation of a degrees the force around x axis 
step_angle=1;
a=0;F_rot=[0;0;0];MMS=[0;0;0];
for i=1:360/step_angle+1
    alpha=step_angle*(i-1)*pi/180;
    a(i)=alpha*180/pi;
    R=[ 1    0            0;...
        0    cos(alpha)  -sin(alpha);...
        0    sin(alpha)   cos(alpha)];
F_rot(:,i)=R*F;
MMS(:,i)=cross(F_rot(:,i),D)*1e3*Time;
end


figure()
set(gcf,'color','w');
plot(a,F_rot(2,:))
hold on
plot(a,F_rot(3,:))
%plot(a,F_rot(1,:))
grid on

MMSx=MMS(1,:);
MMSy=MMS(2,:);
MMSz=MMS(3,:);

figure()
set(gcf,'color','w');
    subplot(3,1,1)
    plot(a,MMSx);
    grid on
    subplot(3,1,2)
    plot(a,MMSy);
    grid on
    subplot(3,1,3)
    plot(a,MMSz);
    grid on

    
    %% Calc maximum momentum 3 new for SCA
clear

d=pi/180;           %deviation in rads of the thruster wrt to x axis, this is the initial in plane xy
Time= 50 * 60;      %time of thrust in seconds (50 minutes)

%the deviation is of d=1 degree
F=0.75*[cos(d);sin(d);0]*1e-3;  %F in [N] 
x_cog=0;        %Distance of CoG from the geometric center x axis in mm
y_cog=3.81;     %Distance of CoG from the geometric center y axis in mm
z_cog=8.23;     %Distance of CoG from the geometric center z axis in mm
x_th=171;       %Position x of point of thrust
y_th=0;         %Position y of point of thrust
z_th=0;         %Position z of point of thrust

%D=[171;y;z]/1e3;    %Distance from the geometric center in m
D=[x_th-x_cog;y_th-y_cog;z_th-z_cog]/1e3;    %Distance from the geometric center in m

%the rotation of a degrees the force around x axis 
step_angle=1;
a=0;F_rot=[0;0;0];MMS=[0;0;0];      %initialize variables
for i=1:360/step_angle+1
    alpha=step_angle*(i-1)*pi/180;
    a(i)=alpha*180/pi;
    R=[ 1    0            0;...
        0    cos(alpha)  -sin(alpha);...
        0    sin(alpha)   cos(alpha)];
F_rot(:,i)=R*F;
MMS(:,i)=cross(F_rot(:,i),D)*1e3*Time;
end

%plot the force as a function of the rotation of the angle a around x axis
figure()
set(gcf,'color','w');
plot(a,F_rot(2,:))
hold on
plot(a,F_rot(3,:))
%plot(a,F_rot(1,:))
grid on


%now plot momentum storage
MMSx=MMS(1,:);
MMSy=MMS(2,:);
MMSz=MMS(3,:);

figure()
set(gcf,'color','w');
    title('Momentum storage in each axis')    
    subplot(3,1,1)
    plot(a,MMSx);
    ylabel('Momentum storage in X')
    xlabel('Angle rotation')
    grid on
    subplot(3,1,2)
    plot(a,MMSy);
    ylabel('Momentum storage in Y')
    xlabel('Angle rotation')
    grid on
    subplot(3,1,3)
    plot(a,MMSz);
    ylabel('Momentum storage in Z')
    xlabel('Angle rotation')
    grid on
    
%% distance of CoG to Thruster point (not useful) just +-1mm in CoG
clear

x_cog=0;        %Distance of CoG from the geometric center x axis in mm
%y_cog=3.81;     %Distance of CoG from the geometric center y axis in mm
y_cog=5;     %Distance of CoG from the geometric center y axis in mm
%z_cog=8.23;     %Distance of CoG from the geometric center z axis in mm
z_cog=5;     %Distance of CoG from the geometric center z axis in mm
x_th=0;       %Position x of point of thrust
y_th=-1;         %Position y of point of thrust
z_th=-1;         %Position z of point of thrust
steps=100;       %in milimeter
sep=2/(steps-1);

for i=1:steps
    for j=1:steps
        D=[x_th-x_cog;y_th(i)-y_cog;z_th(j)-z_cog];    %Distance from the geometric center in m
        dist(i,j)=norm(D);
        z_th(j+1)=sep+z_th(j);
    end
    y_th(i+1)=sep+y_th(i);
end
% D=[x_th-x_cog;y_th(i+1)-y_cog;z_th(i+1)-z_cog]; 
% dist(i+1,i+1)=norm(D);

figure()
    mesh(y_th(1:steps), z_th(1:steps), dist)
    %set(gca,'XTick',[0:30:90]);
    %set(gca,'YTick',[0:30:180]);
    title('Aerodynamic torque profile')
    ylabel('Y [mm]')
    xlabel('X [mm]')
    zlabel('Distance [mm]')
    colorbar
    
%% SAME AS BEFORE< BUT JUST TO PLOT XI  ξ and ψ 
clear

d=pi/180;           %deviation in rads of the thruster wrt to x axis, this is the initial in plane xy
Time= 50 * 60;      %time of thrust in seconds (50 minutes)

%the deviation is of d=1 degree
F=0.75*[cos(d);sin(d);0]*1e-3;  %F in [N] 
x_cog=0;        %Distance of CoG from the geometric center x axis in mm
% y_cog=3.81;     %Distance of CoG from the geometric center y axis in mm
% z_cog=8.23;     %Distance of CoG from the geometric center z axis in mm
x_th=183;       %Position x of point of thrust
y_th=0;         %Position y of point of thrust
z_th=0;         %Position z of point of thrust

%D=[171;y;z]/1e3;    %Distance from the geometric center in m
%D=[x_th-x_cog;y_th-y_cog;z_th-z_cog]/1e3;    %Distance from the geometric center in m

%the rotation of a degrees the force around x axis 
step_angle=1;
a=0;F_rot=[0;0;0];
%MMS=[0;0;0];      %initialize variables
r=8.5;            %radius of circle around YZ plane where the CoG lies

for j=1:360/step_angle+1
    beta=step_angle*(j-1)*pi/180;
    b(j)=beta*180/pi;
    y_cog=r*cos(beta);
    z_cog=r*sin(beta);
    D=[x_th-x_cog;y_th-y_cog;z_th-z_cog]/1e3;    %Distance from the geometric center in m
    
    for i=1:360/step_angle+1    %rotate F around x axis
        alpha=step_angle*(i-1)*pi/180;
        a(i)=alpha*180/pi;
        R=[ 1    0            0;...
            0    cos(alpha)  -sin(alpha);...
            0    sin(alpha)   cos(alpha)];
        F_rot(:,i)=R*F;
        %MMS(i,j,:)=cross(F_rot(:,i),D)*1e3*Time;
        Mom_vec=cross(F_rot(:,i),D)*1e3*Time;   %momentum storage in mNms
        Mom_norm=norm(Mom_vec);
        MMS(i,j,1:3)=Mom_vec;
        MMS(i,j,4)=Mom_norm;
    end
end

%plot the force as a function of the rotation of the angle a around x axis
% figure()
% set(gcf,'color','w');
% plot(a,F_rot(2,:))
% hold on
% plot(a,F_rot(3,:))
% %plot(a,F_rot(1,:))
% grid on


%now plot momentum storage
MMSx=MMS(:,:,1);
MMSy=MMS(:,:,2);
MMSz=MMS(:,:,3);
MMSnorm=MMS(:,:,4);

figure(1)
set(gcf,'color','w');
    subplot(2,2,1)
        mesh(a,b,MMSx);
        ylabel('rotation ξ CoG in XY')
        xlabel('rotation ψ dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in X axis [mNms]') 
        colorbar
    subplot(2,2,2)
        mesh(a,b,MMSy);
        ylabel('rotation ξ CoG in XY')
        xlabel('rotation ψ dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in Y axis [mNms]')   
        colorbar
    subplot(2,2,3)
        mesh(a,b,MMSz);
        ylabel('rotation ξ CoG in XY')
        xlabel('rotation ψ dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in Z axis [mNms]')   
        colorbar
    subplot(2,2,4)
        mesh(a,b,MMSnorm);
        ylabel('rotation ξ CoG in XY')
        xlabel('rotation ψ around X starting in XY')
        grid on
        title('Required Momentum storage norm [mNms]') 
        colorbar
        
% max (MMSy,[],'all')        
% max (MMSz,[],'all')
% max (MMSnorm,[],'all')
%we know that MMS for each wheel is 16.2 mNms,for RWA is 
%1.79*16.2   % =28.998 mNms in each axis 
%taking a 10% margin
%Max mom storage in each axis = 26.3618 mNms
%this is achieved with r<7.5mm, rounding to r=7
%then substracting the 1mm of error in point of thrust:
%r=6mm 
%this is the maximum distance of the CoG from the Geometric Center. in Y-Z

%% now varying two angles
clear

d=pi/180;           %deviation in rads of the thruster wrt to x axis, this is the initial in plane xy
Time= 50 * 60;      %time of thrust in seconds (50 minutes)

%the deviation is of d=1 degree
F=0.75*[cos(d);sin(d);0]*1e-3;  %F in [N] 
x_cog=0;        %Distance of CoG from the geometric center x axis in mm
% y_cog=3.81;     %Distance of CoG from the geometric center y axis in mm
% z_cog=8.23;     %Distance of CoG from the geometric center z axis in mm
x_th=183;       %Position x of point of thrust
y_th=0;         %Position y of point of thrust
z_th=0;         %Position z of point of thrust

%D=[171;y;z]/1e3;    %Distance from the geometric center in m
%D=[x_th-x_cog;y_th-y_cog;z_th-z_cog]/1e3;    %Distance from the geometric center in m

%the rotation of a degrees the force around x axis 
step_angle=1;
a=0;F_rot=[0;0;0];
%MMS=[0;0;0];      %initialize variables
r=8.5;            %radius of circle around YZ plane where the CoG lies

for j=1:360/step_angle+1
    beta=step_angle*(j-1)*pi/180;
    b(j)=beta*180/pi;
    y_cog=r*cos(beta);
    z_cog=r*sin(beta);
    D=[x_th-x_cog;y_th-y_cog;z_th-z_cog]/1e3;    %Distance from the geometric center in m
    
    for i=1:360/step_angle+1    %rotate F around x axis
        alpha=step_angle*(i-1)*pi/180;
        a(i)=alpha*180/pi;
        R=[ 1    0            0;...
            0    cos(alpha)  -sin(alpha);...
            0    sin(alpha)   cos(alpha)];
        F_rot(:,i)=R*F;
        %MMS(i,j,:)=cross(F_rot(:,i),D)*1e3*Time;
        Mom_vec=cross(F_rot(:,i),D)*1e3*Time;   %momentum storage in mNms
        Mom_norm=norm(Mom_vec);
        MMS(i,j,1:3)=Mom_vec;
        MMS(i,j,4)=Mom_norm;
    end
end

%plot the force as a function of the rotation of the angle a around x axis
% figure()
% set(gcf,'color','w');
% plot(a,F_rot(2,:))
% hold on
% plot(a,F_rot(3,:))
% %plot(a,F_rot(1,:))
% grid on


%now plot momentum storage
MMSx=MMS(:,:,1);
MMSy=MMS(:,:,2);
MMSz=MMS(:,:,3);
MMSnorm=MMS(:,:,4);

figure(1)
set(gcf,'color','w');
    subplot(2,2,1)
        mesh(a,b,MMSx);
        ylabel('rotation beta CoG in XY')
        xlabel('rotation alpha dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in X axis [mNms]') 
        colorbar
    subplot(2,2,2)
        mesh(a,b,MMSy);
        ylabel('rotation beta CoG in XY')
        xlabel('rotation alpha dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in Y axis [mNms]')   
        colorbar
    subplot(2,2,3)
        mesh(a,b,MMSz);
        ylabel('rotation beta CoG in XY')
        xlabel('rotation alpha dev. F around X starting in XY')
        grid on
        title('Required Momentum storage in Z axis [mNms]')   
        colorbar
    subplot(2,2,4)
        mesh(a,b,MMSnorm);
        ylabel('rotation beta CoG in XY')
        xlabel('rotation alpha dev. F around X starting in XY')
        grid on
        title('Required Momentum storage norm [mNms]') 
        colorbar
        
% max (MMSy,[],'all')        
% max (MMSz,[],'all')
% max (MMSnorm,[],'all')
%we know that MMS for each wheel is 16.2 mNms,for RWA is 
%1.79*16.2   % =28.998 mNms in each axis 
%taking a 10% margin
%Max mom storage in each axis = 26.3618 mNms
%this is achieved with r<7.5mm, rounding to r=7
%then substracting the 1mm of error in point of thrust:
%r=6mm 
%this is the maximum distance of the CoG from the Geometric Center. in Y-Z