function [y]=MEKFquatf(u)
% Multiplicative Extented discrete-time Kalman filter for attitude 
% estimation based on quaternions and bias rate random walk in the gyro
%
% Global Variables:
%    W6x6:  6x6 covariance matrix of the discrete-time state noise,
%    V:  Matrix of the measurement noise, not used in this implementation,
%    Vmtm:  Covariance of the MTM measurement noise,
%    Vfss:  Covariance of the FSS measurement noise,
%    Vcss:  Covariance of the CSS measurement noise,
%    dt:  sampling time,
%    w_hat: model frequency initial value
%    q_hat0: model quaternion initial value
%    b_hat0: model bias initial value
% Persistent Variables:
%    Xhat: current state estimate = [quaternion;bias]
%    Xhat(1): q0
%    Xhat(2): q1
%    Xhat(3): q2
%    Xhat(4): q3
%    Xhat(5): bias x
%    Xhat(6): bias y
%    Xhat(7): bias z
%    P: current estimation error convariance matrix [6x6]
%    w_hat: angular speed velocity = [wx ; wy ; wz]
%    Xcor: Corrected state variables [7x1]
%    Pcor: Corrected covariance matrix [6x6]
%
% Syntax:
% y=EKFquat(u)
% Input argument:
%   1 * u(1:3):      measurement Gyro
%   2 * u(4:6):      measurement MTM
%   3 * u(7:9):      measurement SS fine
%   4 * u(10:12):    measurement SS coarse
%   5 * u(13):       SS fine enabled
%   6 * u(14):       SS coarse enabled
%   7 * u(15:17):    calculation of magnetic field at position XYZ ECI
%   8 * u(18:20):    calculation of sun direction at position XYZ ECI
%   7 * u(21):       clock with rate of sampling time T (gyro)
% Output argument
%    * y=[Xcor;w_hat;sqrt(diag(Pcor))];
%    * Quaternion (4x1), bias (3x1), angular velocity (3x1), diag cov. (7)

global  W6x6 Vmtm Vfss Vcss dt q_hat0 b_hat0
persistent Xhat P w_hat Xcor Pcor

if u(21)==0  % Initialization on the first measurement, t=0
    Xhat=[1;0;0;0;0;0;0];   %if Xhat is not initialized
    
    Xhat(1)=q_hat0(1); %q0 ... from measurement or TRIAD could be as well?
    Xhat(2)=q_hat0(2);
    Xhat(3)=q_hat0(3);
    Xhat(4)=q_hat0(4);
    Xhat(5)=b_hat0(1); 
    Xhat(6)=b_hat0(2);
    Xhat(7)=b_hat0(3);
    
    w_hat(1)=u(1)-Xhat(5);
    w_hat(2)=u(2)-Xhat(6);
    w_hat(3)=u(3)-Xhat(7);
    w_hat=[w_hat(1);w_hat(2);w_hat(3)];
    
    P=1*eye(6);
    P(1,1)=0.1^2 ; P(2,2)=0.1^2 ; P(3,3) = 0.1^2 ;
    P(4,4)=0.0003^2 ; P(5,5)=0.0003^2 ; P(6,6)=0.0003^2 ;
    
    Xcor=Xhat;
    Pcor=P;
end

if u(21)>0 % when t is positive.
    if u(14) == 1 %at least 2 sensors MTM + CSS
        if u(13) == 1 %3 sensors active: MTM + CSS + FSS
            % to Calculate Kalman Gain:
            % K = P(k|k-1)*H' * (H*P(k|k-1)*H' + R)
            % with H = dh(X)/dX
            b_i=u(15:17);
            s_i=u(18:20);
            q0=Xhat(1);     q1=Xhat(2);     q2=Xhat(3);     q3=Xhat(4);
    
            b_b = (quatrotate([q0 q1 q2 q3], b_i'))';
            s_b = (quatrotate([q0 q1 q2 q3], s_i'))';
    
            H = [2*skew(b_b) , zeros(3);...
                 2*skew(s_b) , zeros(3);...
                 2*skew(s_b) , zeros(3)];  %3 sensors
    
            V =[Vmtm*eye(3),zeros(3,6);zeros(3),Vfss*eye(3),zeros(3);...
                zeros(3,6),Vcss*eye(3)]; 
            % 1. Kalman Gain calculation
            K = P*H' * inv(H*P*H' + V);
            
            % 2. Innovation: [Y(k)-h(X(k))]
            % Correction of state:
            % delta X = K*[Y(k)-h(X(k))] 
            % with Y(k) is the measurement and h(X(k)) the prediction of meas.
            % if K = 0 we only trust the prediction 
            inov =  [u(4:6)-b_b; u(7:9)-s_b; u(10:12)-s_b];
            deltaX = K*inov;
            dq = deltaX(1:3) ; db = deltaX(4:6);
            q_update = quatmultiply( Xhat(1:4)' , [sqrt(1-dq'*dq),dq']  ); 
            beta_update = db + Xhat(5:7);
            Xcor = [q_update';beta_update];
    
            % 3. Correction of Covariance
            % P(k) = (I-K*H) * P(k|k-1) 
            Pcor = (eye(6)-(K*H)) * P ;
        
    else
        % IF there are 2 sensors: MTM + CSS
        % to Calculate Kalman Gain:
        % K = P(k|k-1)*H' * (H*P(k|k-1)*H' + R)
        % with H = dh(X)/dX
        b_i=u(15:17);
        s_i=u(18:20);
        q0=Xhat(1);     q1=Xhat(2);     q2=Xhat(3);     q3=Xhat(4);
    
        b_b = (quatrotate([q0 q1 q2 q3], b_i'))';
        s_b = (quatrotate([q0 q1 q2 q3], s_i'))';
    
        H = [2*skew(b_b) , zeros(3);...
             2*skew(s_b) , zeros(3)];  %2 sensors
        
        V=[Vmtm*eye(3),zeros(3);zeros(3),Vcss*eye(3)]; 
        % 1. Kalman Gain calculation
        K = P*H' * inv(H*P*H' + V);
        
        % 2. Innovation: [Y(k)-h(X(k))]
        % Correction of state:
        % delta X = K*[Y(k)-h(X(k))] 
        % with Y(k) is the measurement and h(X(k)) the prediction of meas.
        inov =  [u(4:6)-b_b; u(10:12)-s_b];
        deltaX = K*inov;
        dq = deltaX(1:3) ; db = deltaX(4:6);
        q_update = quatmultiply( Xhat(1:4)' , [sqrt(1-dq'*dq),dq']  ); 
        beta_update = db + Xhat(5:7); 
        Xcor = [q_update';beta_update];
    
        % 3. Correction of Covariance
        % P(k) = (I-K*H) * P(k|k-1) 
        Pcor = (eye(6)-(K*H)) * P ;
        end
    end
    
    if u(14) == 0   %Eclipse, therefore only MTM available
        % to Calculate Kalman Gain:
        % K = P(k|k-1)*H' * (H*P(k|k-1)*H' + R)
        % with H = dh(X)/dX
        b_i=u(15:17);
        q0=Xhat(1);     q1=Xhat(2);     q2=Xhat(3);     q3=Xhat(4);
    
        b_b = (quatrotate([q0 q1 q2 q3], b_i'))';
           
        H = [2*skew(b_b) , zeros(3)];   %1 sensor
        V=Vmtm*eye(3);
        
        % 1. Kalman Gain calculation
        K = P*H' * inv(H*P*H' + V);
    
        % 2. Innovation: [Y(k)-h(X(k))]
        % Correction of state:
        % delta X = K*[Y(k)-h(X(k))] 
        % with Y(k) is the measurement and h(X(k)) the prediction of meas.
        inov =  [u(4:6)-b_b];
        deltaX = K*inov;
        dq = deltaX(1:3) ; db = deltaX(4:6);
        q_update = quatmultiply( Xhat(1:4)' , [sqrt(1-dq'*dq),dq']  ); %OK?
        beta_update = db + Xhat(5:7); %seems OK?
        Xcor = [q_update';beta_update];
        
        % 3. Correction of Covariance
        % P(k) = (I-K*H) * P(k|k-1) 
        Pcor = (eye(6)-(K*H)) * P ;
    end
    
    %NOW PREDICTION
     % 4. Correction of angular velocity
    % estimated angular velocity = measurement - estimated bias
    w_hat = u(1:3)-Xcor(5:7);
    
    % 5. State prediction:
    % q- k+1 = THETA(w)* q+ k
    q=[Xcor(2:4);Xcor(1)];
    normW = norm(w_hat);
    psi=(sin(0.5*normW*dt)*w_hat)/normW;
    
    THETA = [ cos(0.5*normW*dt)*eye(3)-skew(psi) , psi ;...
                    -psi'      ,      cos(0.5*normW*dt)];      
    
    q=THETA*q;
    Xhat=[q(4);q(1);q(2);q(3);Xcor(5:7)];   %Prediction DONE!
      
    % Computation on the linearized dynamics Ad or PHI
    PHI11 = eye(3)-skew(w_hat)*(sin(normW*dt))/normW+...
            ((skew(w_hat))^2)*(1-cos(normW*dt))/normW^2;
    PHI12 = skew(w_hat)*(1-cos(normW*dt))/normW^2 - eye(3)*dt-...
            (skew(w_hat))^2*(normW*dt-sin(normW*dt))/normW^3;
    PHI21=zeros(3);
    PHI22=eye(3);
    PHI = [ PHI11 , PHI12 ; PHI21 , PHI22];
             
    % 6. Prediction error covariance matrix:
    % P(k|k-1) = Φ * P(k-1) * Φ' + G*Q*G'  (GQG' is W6x6)
    P = PHI*Pcor*PHI' + W6x6;
    
end
% Output:
 y=[Xcor;w_hat;sqrt(diag(Pcor))];
end