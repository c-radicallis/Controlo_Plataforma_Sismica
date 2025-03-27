%% Structure parameters

mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;

% 1st mode
m1 = mass; % kg
f1 = 0.4; % Hz   % 1.5 < f1 < 4
zeta1 = 0.1 ; % 2 < zeta1 < 10
%2nd mode
m2 = mass; % kg
f2 =3; % Hz % 6 < f2 < 10
zeta2 = 0.06; % 5 < zeta2 < 25r

% Controller
k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);

%coupled 2DOF system
c1 = zeta1*2*m1*2*pi*f1; %N/m/s
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
digits(1e5)
AA = vpa([-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
          k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
              0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
              0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
              0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2]);
          
BB = vpa([k_svk_q/tau_sv ; zeros(7,1)]);
CC = vpa([zeros(1,2), 1 , zeros(1,5)]);  % measuring xT
DD = 0;

%%
clear
%'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\
load('mat and fig files\obsv_ctrb_vpa1e5.mat')

% ss_model =
% 
%   A = 
%                x1          x2          x3          x4          x5          x6          x7          x8
%    x1      -40.65           0           0           0           0           0           0           0
%    x2    3.63e+12  -4.878e+04           0           0           0  -4.521e+10           0           0
%    x3           0           0           0           0           0           1           0           0
%    x4           0           0           0           0           0           0           1           0
%    x5           0           0           0           0           0           0           0           1
%    x6           0   0.0005063      -13.03       13.03           0      -3.435       0.509           0
%    x7           0           0       12.87      -187.2       174.4      0.5027      -2.765       2.262
%    x8           0           0           0       174.4      -174.4           0       2.262      -2.262
% 
%   B = 
%             u1
%    x1  0.07864
%    x2        0
%    x3        0
%    x4        0
%    x5        0
%    x6        0
%    x7        0
%    x8        0
% 
%   C = 
%        x1  x2  x3  x4  x5  x6  x7  x8
%    y1   0   0   1   0   0   0   0   0
% 
%   D = 
%        u1
%    y1   0


%%
obs = vpa(obsv(AA, CC));
r_obsv = rank(obs)
obs=double(obs)
ctrlb = vpa(ctrb(AA,BB));
r_ctrlb = rank(ctrlb)
ctrlb = double(ctrlb)

% r_obsv =
% 
%      8
% 
% 
% obs =
% 
%    1.0e+28 *
% 
%          0         0    0.0000         0         0         0         0         0
%          0         0         0         0         0    0.0000         0         0
%          0    0.0000   -0.0000    0.0000         0   -0.0000    0.0000         0
%     0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000
%    -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000
%     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000
%    -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000
%     1.0118   -0.0000   -0.0000    0.0000   -0.0000   -0.0126    0.0000   -0.0000
% 
% 
% r_ctrlb =
% 
%      8
% 
%  ctrlb =
% 
% 1.0e+39 *
% 
% 0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000
%      0    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0001    3.6671
%      0         0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000
%      0         0         0         0    0.0000   -0.0000    0.0000   -0.0000
%      0         0         0         0         0    0.0000   -0.0000    0.0000
%      0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000
%      0         0         0    0.0000   -0.0000    0.0000   -0.0000    0.0000
%      0         0         0         0    0.0000   -0.0000    0.0000   -0.0000


%%
clear
%'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\
load('mat and fig files\obsv_ctrb_vpa1e5.mat')

nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

% Define process noise matrices
G = eye(nx);         % Process noise matrix (assuming full-state noise)
H = zeros(ny, nx);   % No direct noise feedthrough

% Create the augmented system
digits(1e4)
sys_aug = ss(double(AA), [double(BB) G], double(CC), [double(DD) H]);

% Assign input groups:
% - Channel 1 is control,
% - Channels 2 to (nx+1) are noise.
sys_aug.InputGroup.control = 1;
sys_aug.InputGroup.noise   = 2:(nx+1);
% Explicitly set the KnownInput group to only the control input
sys_aug.InputGroup.KnownInput = 1;

% Design the LQI controller for the original system
Q = blkdiag(eye(nx), eye(ny));
R = eye(nu);
K = lqi(ss(double(AA), double(BB),  double(CC),  double(DD)), Q, R)

% Define noise covariance data
% Here Qn should be for process noise (nx-by-nx) and Rn for measurement noise (ny-by-ny)
Qn = eye(nx);
Rn = eye(ny);

% Construct the Kalman estimator using the augmented system
kest = kalman(sys_aug, Qn, Rn)

% Now, lqgtrack expects that the number of rows in K (control actions) matches 
% the number of known input channels in kest (which is now 1).
trksys = lqgtrack(kest, K)

