clear all

%Servo-valve parameters
tsv=0.0246      %Valve time constant (tsv=0.0246 s)
ksvkq=1934.5    %Valve flow gain (ksvkq=1934.5 cm3/s/V)
kpl=1.67401e-7  %Valve pressu
% re gain & leakadge factor (kpl=1.67401e-7 m3/s/kPa)
Be=193716.28    %Oil Bulk modulus (Be=193716.28 kPa)
Vt=0.002659     %Oil Volume on actuator chamber (Vt=0.002659 m3)
A=0.012456      %Piston area (A=0.012456 m2)
mp=1.9751       %Platen mass (mp=1.9751 t)
ctm=5.78        %Total damping, actuator + platen (ct=5.78 kN s/m1)
kh=4*Be*A^2/Vt %(kPa m1)

%Additional 1DOF model
msp=4                             %Rigid mass (t)
fn1dof=2                          %Natural frequency (Hz)
zeta1dof=0.02                     %damping ratio (2%)
wn1dof=2*pi*fn1dof           %Natural frequency (rad/s)
cct=2*msp*wn1dof                  %critical damping (t/s)
csp=zeta1dof*cct                  %damping (t/s)
ksp=wn1dof^2*msp                 %stiffness (t/s^2)


%Controller
kp=1.2993       %Pgain (kp=1.2993 V/cm)

s=tf('s')

%1DOF model
Hsp=-msp*s*s/(msp*s*s+csp*s+ksp)
mt=mp+msp*(1+Hsp)
%mt=mp+msp      %massa rigida


Gxr=10^-4*A*kh*kp*ksvkq/(10^-4*A*kh*kp*ksvkq + s*(s*tsv + 1)*(A^2*ctm*s + A^2*kh + A^2*mt*s^2 + kpl*ctm*kh + kpl*kh*mt*s))

%FT do anel aberto Gaa=xp/isv
Gaa=10^-4*A*kh*ksvkq/(s*(s*tsv + 1)*(A^2*ctm*s + A^2*kh + A^2*mt*s^2 + kpl*ctm*kh + kpl*kh*mt*s))

%FT do anel fechado
Gxraf=Gaa*kp/(1+Gaa*kp)


%%%  2DoF system
% System
m2   = 100;       %      [kg]
m1   = 2000;         %       [kg]
k2  = 5000;       %        [N/m]
k1  = 200000;     %        [N/m]
c2  = 4000;          %      [N.s/m]
c1=2000;   %      [N.s/m]



%  State space model

% M=[ m1 0;
%         0    m2];
% 
% C=[ c1+c2  -c2;
%             -c2    c2];
% 
% K=[k1+k2  -k2;
%         -k2       k2];

n=2;
M=[ m1 m2];
C=[ c1  c2];
K=[k1 k2];

[A,B]= Model_nDOF_Spring_Mass_Damper_SIxOsystem_Force(n,M,C,K);

% A =
%   0	0	1	0
%   0	0	0	1
%   -100	100	-1	1
%   2000	-2050	20	-60

%  B = 
% 0
% 0
% 0.0005
% 0

%this does not consider the damper c1
% A = [ 0                        1              0               0       ;
%       -(k2+k1)/m1      -c2/m1     k2/m1    c2/m1    ;
%        0                           0               0              1       ;
%       k2/m2            c2/m2      -k2/m2   -c2/m2   ];
% 
% B = [ 0     ;
%       k1/m1  ;
%       0     ;
%       0     ];

C = [ 1 0 0 0 ; 
         0 0 1 0 ];
D = [0 ; 0];

sys = ss(A,B,C,D);


%%

dados = load('elcentro.txt');

t_vector = dados(:,1);
t_step = t_vector(2)
ddx = [t_vector dados(:,2)];
ddy = [t_vector  dados(:,3)];

%%

% EspetroResposta_modded %( 'elcentro')