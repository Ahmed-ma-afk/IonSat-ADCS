%% Reaction wheel configurations
clear
clc
close all

%cases = cell(1,2); %Case 1 is the example, case 2 is the IonSat RWA
cases = cell(1); %Case 1 is the IonSat RWA

% cases{1}.A = [...
%     1, -1, 1, -1;
%     1, 1, -1, -1;
%     1, 1, 1, 1];
% cases{1}.Tmax = [100,100,100,100];


% alpha = 26.6 * pi/180;      %from deg to rad
% sat.wheel.repartition_matrix_4RW = [[cos(alpha),      0,    -cos(alpha),      0    ];...
%                                     [     0,     cos(alpha),     0,     -cos(alpha)];...
%                                     [sin(alpha), sin(alpha), sin(alpha), sin(alpha)]];
                                
alpha = 26.6 * pi/180;      %from deg to rad
cases{1}.A = [[cos(alpha),      0,    -cos(alpha),      0    ];...
             [     0,     cos(alpha),     0,     -cos(alpha)];...
             [sin(alpha), sin(alpha), sin(alpha), sin(alpha)]];
cases{1}.Tmax = [1,1,1,1];

% cases{3}.A = [...
%     1, -1, 1, -1, 1, 0, 1, 1;
%     1, 1, -1, -1, 0, 1, 1, -1;
%     1, 1, 1, 1, 0, 0, 0, 0];
% cases{3}.Tmax = [10,20,30,40,50,60,70,80];

for i = 1:1:length(cases)
    A = cases{i}.A;
    RWN = size(A, 2);
    for j = 1:1:RWN
        A(:,j) = A(:,j) / norm(A(:,j));
    end
    W = A * diag(cases{i}.Tmax);
    cases{i}.W = W;
end

%% Draw the torque envelopes
for c = 1:1:length(cases)
    W = cases{c}.W;
    RWN = size(W, 2);
    [C, K] = CplnrWhlsCombn(W);

    [X,Y,Z] = sphere(128);
    N = size(X,1);

    for i = 1:1:N
        for j = 1:1:N
            v = [X(i,j);Y(i,j);Z(i,j)];
            [r, Tmax_L2] = MaxRwTorque(v, W, C);
            v = W * r;
            X(i,j) = v(1);
            Y(i,j) = v(2);
            Z(i,j) = v(3);
        end
    end
    
    figure()
    set(gcf,'color','w');
    surf(X,Y,Z, 'EdgeColor', 'k', 'FaceColor', 'w')
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title(['Torque Envelope [mNm] for ', num2str(RWN), ' RWs'])

end

%Now to plot in Spherical Coordinates
% for m=1:(128+1)
%     for n=1:(128+1)
%         [azimuth(m,n),elevation(m,n),Tmax(m,n)] = cart2sph(X(m,n),Y(m,n),Z(m,n));
%     end
% end
[azimuth,elevation,Tmax] = cart2sph(X,Y,Z);
figure()
set(gcf,'color','w');
%surf(azimuth*180/pi,elevation*180/pi,Tmax, 'EdgeColor', 'k', 'FaceColor', 'w')
%mesh(azimuth*180/pi,elevation*180/pi,Tmax)
surf(azimuth*180/pi,elevation*180/pi,Tmax)
xlabel('Azimuth')
ylabel('Elevation')
zlabel('Torque max. [mNm]')
title(['Torque Envelope for IonSat RWA'])
colorbar