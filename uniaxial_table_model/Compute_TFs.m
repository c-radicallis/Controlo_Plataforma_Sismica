function [s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT,G_Fp_isv ,c1,c2,k1,k2 , AA , BB , CC , DD ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2)
% Compute_TFs  Frequency-domain and state-space Transfer Functions for a
%               coupled 2-DOF system with a hydraulic servo-valve.
%
% SYNTAX:
%   [s, G_T, G_1, G_2, G_T1, G_21, G_svq, G_csv, ...
%    G_x2_x1, G_x1_xT, G_xT_Fp, G_Fp_xref, G_xT_xref, ...
%    G_x1_xref, G_x2_xT, G_Fp_isv, c1, c2, k1, k2, AA, BB, CC, DD] = ...
%    Compute_TFs(G_c, mT, cT, m1, m2, f1, zeta1, f2, zeta2)
%
% INPUTS:
%   G_c    — controller transfer function (TF) from valve to current (–)
%   mT     — platen mass [kg]
%   cT     — total damping (actuator + platen) [N·s/m]
%   m1     — mass of first DOF [kg]
%   m2     — mass of second DOF [kg]
%   f1     — natural frequency of mode-1 [Hz]
%   zeta1  — damping ratio of mode-1 (–)
%   f2     — natural frequency of mode-2 [Hz]
%   zeta2  — damping ratio of mode-2 (–)
%
% OUTPUTS:
%   s          — Laplace operator (for defining TFs)
%   G_T        — TF of the platen (mT*s^2 + (cT+c1)*s + k1)
%   G_1        — TF of the coupled 1st mass (m1*s^2 + (c1+c2)*s + k1+k2)
%   G_2        — TF of the coupled 2nd mass (m2*s^2 + c2*s + k2)
%   G_T1       — coupling TF between platen and mass-1 (c1*s + k1)
%   G_21       — coupling TF between mass-1 and mass-2 (c2*s + k2)
%   G_svq      — servo-valve TF (k_svk_q/(1 + τ_sv·s))
%   G_csv      — combined controller+servo-valve TF (G_c·G_svq)
%   G_x2_x1    — displacement TF: x2/x1
%   G_x1_xT    — displacement TF: x1/xT
%   G_xT_Fp    — transfer from piston flow Fp to platen xT
%   G_Fp_xref  — flow-to-reference TF including leakage & bulk modulus
%   G_xT_xref  — closed-loop TF: xT to reference
%   G_x1_xref  — closed-loop TF: x1 to reference
%   G_x2_xT    — closed-loop TF: x2 to xT
%   G_Fp_isv   — TF: isv (servo-valve current) to piston flow Fp
%   c1, c2     — damping coefficients for the two modes [N·s/m]
%   k1, k2     — stiffness coefficients for the two modes [N/m]
%   AA, BB, CC, DD
%              — state-space matrices of the original plant model
%
% EXAMPLE:
%   % Define parameters
%   Gc    = tf(1, [0.01 1]);
%   [s,GT,~,~,~,~,~,~,~,~,~,~,~,~,~,~,c1,c2,k1,k2,AA,BB,CC,DD] = ...
%     Compute_TFs(Gc, 1975.1, 5780, 500, 300, 5, 0.02, 15, 0.03);
%
%   % Plot the closed-loop TF from reference to platen xT
%   T = feedback(G_Fp_xref*G_xT_Fp, 1);
%   bode(T);
%
% See also TF, SS, vpasolve, tf, ss.

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

%% Original Plant Definition
AA = [-1/tau_sv, 0              , 0      , 0          , 0     , 0              , 0         , 0     ;
          k_h/A , -k_h*k_pl/(A^2), 0      , 0          , 0     , -k_h           , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 1              , 0         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 1         , 0     ;
              0 ,              0 , 0      , 0          , 0     , 0              , 0         , 1     ;
              0 ,           1/mT , -k1/mT , k1/mT      ,  0    , (-cT - c1) /mT ,  c1 /mT   , 0     ;
              0 ,            0   , k1/m1  ,(-k1-k2)/m1 , k2/m1 , c1/m1          ,(-c1-c2)/m1, c2/m1 ;
              0 ,            0   , 0      , k2/m2      , -k2/m2, 0              , c2/m2     , -c2/m2];        
BB = [k_svk_q/tau_sv ; zeros(7,1)];
CC = [zeros(1,2), 1 , zeros(1,5)];  % measuring xT
DD = 0;

end