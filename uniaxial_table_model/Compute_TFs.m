function [s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT,G_Fp_isv ,c1,c2,k1,k2 , ss_model ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2)

%coupled 2DOF system
%  structure parameters
% 1st mode
% m1 = mass; % kg
% f1 = 2; % Hz
% zeta1 = 0.02 ; 

% % k1 = m1*(2*pi*f1)^2; %N/m
c1 = zeta1*2*m1*2*pi*f1; %N/m/s

% %2nd mode
% m2= m1; % kg
% f2 = 10; % Hz
% zeta2 = 0.05; %

c2 = zeta2*2*m2*2*pi*f2; %N/m/s

syms k1 k2
assume(k1 ,"positive")
assume(k2 ,"positive")

Y = vpasolve([
    (m1*k2 + m2*(k1+k2))/(2*m1*m2) - 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f1)^2  ,
    (m1*k2 + m2*(k1+k2))/(2*m1*m2)+ 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f2)^2,
], [k1,k2]);

k1 = double(Y.k1);
k2 = double(Y.k2);

% sprintf("k1 = %e",k1)
% sprintf("k2 = %e",k2)

%Servo-valve parameters
% All converted to SI units
tau_sv=0.0246     ; %Valve time constant (tsv=0.0246 s)
k_svk_q=1934.5*(1e-2)^3   ; %Valve flow gain (ksv∙kq=1934.5 cm3/s/V)

%from Gidewon k_pl =  k_c + C_l 
k_pl=1.67401e-7/1e3 ; %Valve pressure gain & leakadge factor (kpl=1.67401e-7 m3/s/kPa)

Be=193716.28*1e3 ;  %Oil Bulk modulus (Be=193716.28 kPa)
Vt=0.002659   ;  %Oil Volume on actuator chamber (Vt=0.002659 m3)
A=0.012456    ;  %Piston area (A=0.012456 m2)
k_h=4*Be*A^2/Vt*1e3; %(kPa m1)

% mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
% cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)


% TF's
s=tf('s')  ;

G_T=mT*s^2+(cT+c1)*s+k1;
G_1=m1*s^2+(c1+c2)*s+k1+k2;
G_2=m2*s^2+c2*s+k2;

G_T1=c1*s+k1;
G_21=c2*s+k2;

G_svq = k_svk_q/(1+tau_sv*s);

G_csv=G_c*G_svq;

G_x2_x1= G_21/G_2;
G_x1_xT = G_T1*G_2/(G_1*G_2-G_21^2);
G_xT_Fp = (G_1*G_2-G_21^2)/(G_T*G_1*G_2-G_T*G_21^2-G_2*G_T1^2);

G_Fp_xref = G_csv/( k_pl/A + A*s/k_h +G_xT_Fp*(G_csv + A*s));

G_xT_xref = G_Fp_xref * G_xT_Fp;
G_x1_xref = G_x1_xT*G_xT_xref;
G_x2_xT = G_x2_x1 * G_x1_xT;

G_Fp_isv = A*G_svq/( k_pl+A^2*s/k_h+A^2*s*G_xT_Fp );

% state space model

% Define the Mass matrix M
% M = [mT, 0,   0;
%      0,   m1, 0;
%      0,   0,   m2];

% Define the Damping matrix C
% C = [cT + c1, -c1,       0;
%      -c1,      c1 + c2, -c2;
%      0,        -c2,       c2];
% 
% Define the Stiffness matrix K
% K = [k1, -k1,  0;
%      -k1, k1 + k2, -k2;
%      0,   -k2,  k2];

%% Original Plant Definition
A = vpa([-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
          k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
              0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
              0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
              0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2],5000);
          
B = vpa([k_svk_q/tau_sv ; zeros(7,1)],5000);
C = eye(8); %vpa([zeros(1,2), 1 , zeros(1,5)],5000);  % measuring xT
D = 0;

ss_model = ss(double(A), double(B), double(C), D);

obs = obsv(double(A), double(C));
r_obsv = rank(obs)
ctrlb = ctrb(double(A), double(B));
r_ctrlb = rank(ctrlb)

%% LQI Design (using original plant)
nx = size(ss_model.A,1);  % 8 states
ny = size(ss_model.C,1);  % 1 output
nu = size(ss_model.B,2);

% Q weighting matrix for the augmented system (8+1 = 9 states)
Q_lqi = blkdiag(eye(nx), eye(ny));  % 9x9 matrix
R = eye(nu);

% lqi will internally augment ss_model (adds an integrator) to form a 9-state system.
K = lqi(ss_model, Q_lqi, R)

%% Kalman Estimator Design for Tracking
% Augment the plant manually to include an integrator:
A_aug = [double(A), zeros(nx,ny);
         -double(C), zeros(ny,ny)];
%B_aug2 = [zeros(nx,ny); eye(ny)];  % extra input channel for the reference
% Form the estimator plant with two inputs:
B_aug = [double(B), zeros(nx,nu);
                   0            ,   1];  
C_aug = [double(C), eye(ny,ny)];
D_aug = zeros(ny, 2);

ss_model_kalman = ss(A_aug, B_aug, C_aug, D_aug);

% Noise covariance matrices (tune as needed)
Qn_aug = eye(size(B_aug,2));  % Process noise covariance (9x9)
Rn_aug = 1;           % Measurement noise covariance

kest = kalman(ss_model_kalman, Qn_aug, Rn_aug)

%% Connect Estimator and State-Feedback Gain to Form LQG Servo Controller
% 'lqgtrack' requires that the estimator (kest) has at least 2 inputs and 
% an appropriate number of outputs (state estimates) that match the augmented plant.
trksys = lqgtrack(kest, K, "2dof")

end