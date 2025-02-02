function [s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT,G_Fp_isv  ]=Compute_TFs(G_c)

%coupled 2DOF system
%  structure parameters
% 1st mode
m1= 1e3; % kg
f1 = 2; % Hz
zeta1 = 0.02 ; 

% k1 = m1*(2*pi*f1)^2; %N/m
c1 = zeta1*2*m1*2*pi*f1; %N/m/s

%2nd mode
m2= 1e3; % kg
f2 = 10; % Hz
zeta2 = 0.05; %

c2 = zeta2*2*m2*2*pi*f2; %N/m/s

syms k1 k2
assume(k1 ,"positive")
assume(k2 ,"positive")

Y = vpasolve([
    (m1*k2 + m2*(k1+k2))/(2*m1*m2) - 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f1)^2  ,
    (m1*k2 + m2*(k1+k2))/(2*m1*m2)+ 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f2)^2,
], [k1,k2])

k1 = double(Y.k1)
k2 = double(Y.k2)


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

mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)


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
G_xT_Fp = G_1*G_2-G_21^2/(G_T*G_1*G_2-G_T*G_21^2-G_2*G_T1^2);

G_Fp_xref = G_csv/( k_pl/A + A*s/k_h +G_xT_Fp*(G_csv + A*s));

G_xT_xref = G_Fp_xref * G_xT_Fp;
G_x1_xref = G_x1_xT*G_xT_xref;
G_x2_xT = G_x2_x1 * G_x1_xT;

G_Fp_isv = A*G_svq/( k_pl+A^2*s/k_h+A^2*s*G_xT_Fp );

end