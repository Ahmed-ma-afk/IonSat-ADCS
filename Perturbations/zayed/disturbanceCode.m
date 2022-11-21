clear
Cdrag=2.6;%%drag_coefficient
norme2=6*1e7;%Velocity square
%gravity center
OGX=2.5*1e-3;
OGY=1*1e-2;
OGZ=1*1e-2;
G=[OGX,OGY,OGZ];
%density_of_the_atmosphere
d=1e-11;

step=1;
%%%spherique_orientation.
theta=(0:step:180)*pi/180;%%the_angle_(z,v)_v_is_the_velocity_vector
phi=(0:step:360)*pi/180;%%the_angle_(x,pv)_pv_is_the_projection_of_v_in_XY_paln.
T=zeros(length(theta),length(phi));


for itheta=1:length(theta)
    itheta_ANG=theta(itheta);
    for iphi=1:length(phi)
        iphi_ANG=phi(iphi);
        tmp=calctorque(itheta_ANG, iphi_ANG,G);
        tmp=-0.5*d*norme2*Cdrag*tmp;%%torque_in_body_frame
        T(itheta,iphi)=sqrt(tmp(1)^2+tmp(2)^2+tmp(3)^2);%%norme_of_the_torque.
    end
end
[Phi,Theta]=meshgrid(phi,theta);
mesh(Phi*180/pi,Theta*180/pi,T)
%mesh(Phi*180/pi,Theta*180/pi,log10(T))
%surf(Phi*180/pi,Theta*180/pi,T)
xlabel('Phi en degré')
ylabel('Theta en degré')
title('Norme du couple aerodynamique')