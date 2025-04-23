clear; close all;

% %% --- original and inverse TF
% G    = tf([1 2], [1 5 6]);
% Ginv = inv(G);
% %--- create input and convolve
% fs  = 100; 
% t_u = 0:1/fs:10;
% u   = sin(2*pi*0.5*t_u);

%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 4; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =10; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller
% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2, AA,BB,CC,DD
[s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~,~,~,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);

lim_displacement = 0.1; % m % Limits
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

s=tf('s') ;
x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref);
scale=1;
while max_xref > lim_displacement % Scaling down if necessary
    scale = 0.95*scale;
    ddx_ref = 0.95*ddx_ref;
    x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
    max_xref = max(x_ref);
end
v_ref =  lsim(1/s,  ddx_ref , t_vector ,'foh');
max_vref = max(v_ref);

% tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');% Use PID tuner app to generate a PID controller for the system
% G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
% [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ ,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

t_u = t_vector';
u   = x_ref';
G = G_xT_xref ;
Ginv = inv(G_xT_xref );

%%
%--- frequency‐grid
fs  = 1/t_step;           % Hz
N   = 2^12;          
f   = (-N/2:N/2-1)*(fs/N);
w   = 2*pi*f;

%--- sample H(jw) and get h[n]
H   = squeeze(freqresp(Ginv, w));
h   = ifft(ifftshift(H));
t_h = (-N/2:N/2-1)/fs;


y_full = conv(u, h, 'full')*(1/fs);

%--- extract causal portion
startSample = find(abs(t_h) < 1e-6,1);            % center of h
y = y_full(startSample : startSample + numel(u)-1);
plot(t_u , u);
hold on
plot(t_u, y);
xlabel('Time (s)'); ylabel('y(t)'); title('Inverse‐system response');

% %"Zero-phase” filtering with filtfilt
% [b,a] = tfdata(Ginv,'v');    % get numerator and denominator
% y = filtfilt(b,a, u);         % forward–backward filter → zero-phase, non-causal
% plot(t_u,y)
% %Error using filtfilt>getCoeffsAndInitialConditions
% %   First denominator filter coefficient must be nonzero.