function [T]=calctorque(theta, phi,G)
%%%%%%%
AX=0.022;%%surface_facing_X_direction
AY=0.034;%surface_facing_Y_direction   %%in m2
AZ=0.0748;%surfacez_facing_Z_direction
OX=0.17;
OY=0.11;
OZ=0.05;
%about deployable_solar_plan
CRS=[0,0.195,0.05];%%the_center_of_the_right_one.
CLS=[0,-0.195,0.05];%%the_center_of_the_left_one.
ESS=0.0578;%%surface_of_one
AZES=0.1156;%surface_of_the_2
%about_the_surface_facing Y
Cmy=[0,-0.11,0];%the_center_of_the_one_facing_-Y.
%%velocity_vector
v=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];%%here_we_use_spherique_oriatation_and_v_represente_the_normalize_vector_of_the_velocity.

u1=v(3);%%is_the_dot_product_between_Z_and_v
u2=v(2);%%is_the_dot_product_between_y_and_v
torque=[0,0,0];
psi=1e-10;%% optimisation

if u1 > psi&&u2>psi    %%u1>0 and u1 !=1  and u2>0
    torque=projectionTorque1(G,v);%%here_we_calculate_the_"torque"_due_to_the_aire_facing_Y_direction_and_the_wind.(caused_by_the_right_deployable_solar_plan)
    Yshadow=u2;
    Zh=u1;
    torque=torque*Yshadow;
    torque=torque+(cross(CRS-G,v)+cross (CLS-G,v))*ESS*Zh;%%+torque_due_to_the_deployable_solar
    
elseif u1<-psi&&u2<-psi %% u1<0 and u1!=-1 and u2<0
    torque=projectionTorque2(G,v);%%here_we_calculate_the_"torque"_due_to_the_aire_of_the_right_deployable_solar_plan_facing_-Z_and_the_wind.(caused_by_the_aire_facing_Y)
    Ysh=u2;
    Zh=u1;
    torque=torque*abs(Zh);
    torque=torque+ cross(CLS-G,v)*ESS*abs(Zh)+cross(Cmy-G,v)*AY*abs(Ysh);%%+torque_due_to_the_left_deployable_solar_plan+torque_due_to_the_surface_facing_-Y.
    
elseif u1>psi&&u2<-psi %% u1>0 and u1!=1 and u2 <0
    torque=projectionTorque3(G,v);%%here_we_calculate_the_"torque"_due_to_the_aire_facing_-Y_direction_and_the_wind.(caused_by_the_left_deployable_solar_plan)
    Ysh=u2;
    Zh=u1;
    torque=torque*abs(Ysh);
    torque=torque+(cross(CRS-G,v)+cross(CLS-G,v))*Zh*ESS;%%+torque_due_to_the_deployable_solar
    
    
elseif u1<-psi&&u2>psi %% u1<0 and u1 !=-1 and u2>0
    torque=projectionTorque4(G,v);%%here_we_calculate_the_"torque"_due_to_the_aire_of_the_left_deployable_solar_plan_facing_-Z_and_the_wind.(caused_by_the_aire_facing_-Y)
    Ysh=u2;
    Zh=u1;
    torque=torque*abs(Zh);
    torque=torque+cross(CRS-G,v)*ESS*abs(Zh)-cross(Cmy+G,v)*AY*Ysh;%%+torque_due_to_the_right_deployable_solar_plan+torque_due_to_the_surface_facing_Y.
     
else
    z=[0,0,1];
    Zh=u1;
    Yh=u2;
    t1=2*(cross(z,v)*OZ-cross(G,v))*abs(Zh)*ESS;%%torque_due_to_the_deployable_solar_plan.
    y=[0,1,0];
    t2=(cross(y,v)*Yh*OY-abs(Yh)*cross(G,v))*AY;%%torque_due_to_the_surface_in_-Y_or_Y_direction.(depending_on_the_sign_of_X_dot_v)
    torque=t1+t2;
end
x=[1,0,0];
z=[0,0,1];
torque=torque+(cross(x,v)*v(1)*OX-abs(v(1))*cross(G,v))*AX;%%torque_due_to_the_surface_facing_X_or_-X_direction.()
torque=torque+(cross(z,v)*v(3)*OZ-abs(v(3))*cross(G,v))*AZ;%%torque_due_to_the_surface_facing_Z_or_-Z_direction.(with_out_deployable_solar_plan).
% test
% y=[0,1,0];
% torque=torque+(cross(y,v)*v(2)*OY-abs(v(2))*cross(G,v))*AY;%%if there
% %ware no solare deployable solar panals.....

T=torque;
end

