clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

%% Loading Model with Standard Tune
mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 1.5; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =6; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

c1 = zeta1*2*m1*2*pi*f1; c2 = zeta2*2*m2*2*pi*f2; %N/m/s%coupled 2DOF system

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2,AA , BB , CC , DD 
[~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , AA , BB , CC , DD  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

Gmini=tf(100*2*pi,[1 100*2*pi]);

fig1 = figure(1); %clf; % Clear figure if needed
ax1 = axes(fig1); hold(ax1, 'on');
opts1=bodeoptions('cstprefs');
opts1.FreqUnits = 'Hz';
opts1.XLim={[1 40]};
% opts1.YLim={[-40 1]};
% opts1.MagVisible='off';

bodeplot(G_xT_xref,opts1);
hold on
bodeplot(Gmini,opts1);
bodeplot(Gmini*G_xT_xref,opts1);
%bodeplot( tf(10*2*pi,[1 10*2*pi])*tf(10*2*pi,[1 10*2*pi])*tf(10*2*pi,[1 10*2*pi]) ,opts1);
legend();
grid on;

