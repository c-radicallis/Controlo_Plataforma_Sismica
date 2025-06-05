clear;clc;close all;

%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

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

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2, ss_model 
[~,~,~,~,~ ,~ ,~,~,~,~,~,~,G_xT_xref,~,~ , ~ ,~,~,~,~ , ~ , ~ , ~ , ~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

%% 
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'\Adapting_Driver_Signal\PRJ_project\
filename = 'pink_noise_40Hz_T3mm.drv.txt'; 
loadTXT(filename)

%%   % --- Simulation --
t_vector = drv_time_vector_pink_noise_40Hz_T3mm;
t_step = t_vector(2);
x_drv    = drvPosA1T_pink_noise_40Hz_T3mm;
ddx_drv = secondDerivativeTime(x_drv,t_step);

x_acq = lsim(G_xT_xref ,  x_drv ,t_vector,'foh');
ddx_acq = secondDerivativeTime(x_acq,t_step);
%%

 hold on; grid on; legend()
 plot(t_vector , x_drv)
 plot(t_vector,x_acq)

%%  --- Write .txt file from simulated "acquired" data 
filename_acq = strrep(filename, '.drv.txt', '_acq.txt');
 writeTXT(t_vector , x_acq , ddx_acq , filename_acq)

% Run .txt to .LTF conversion
% Run .LTF to .txt conversion

%% Load the file converted-deconverted file and compare to the original

filename = 'pink_noise_40Hz_T3mm_acq.txt.acq.txt' ;
loadTXT(filename)
 plot(acq_time_vector_pink_noise_40Hz_T3mm_acq_txt , acqPosT_pink_noise_40Hz_T3mm_acq_txt,'.--')