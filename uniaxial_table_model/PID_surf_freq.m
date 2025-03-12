clear 
close all
%clc

%%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);
ddx = [t_vector dados(:,2)];

% Limits
lim_displacement = 100; % mm
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

% Scaling down if necessary
% displacements in milimeters
s=tf('s');
x_ref = lsim(1e3/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref)
scale=1;
while max_xref > lim_displacement
    scale = 0.95*scale;
    ddx_ref = 0.95*ddx_ref;
    x_ref = lsim(1e3/s^2,  ddx_ref , t_vector ,'foh');
    max_xref = max(x_ref);
end
scale
max_xref
v_ref =  lsim(1/s,  ddx_ref , t_vector ,'foh');
max_vref = max(v_ref)

%% Finding Response Spectre of Ground
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 200;
%f_vector = linspace( f_i , f_n , n_points); 
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

[picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , ddx_ref , f_vector );

%% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;

% Before the loops, initialize arrays to store f1, f2, and mse values
f1_arr = [];
f2_arr = [];
mse_arr   = [];
isv_arr=[];
Fp_arr=[];


%%
freq_list =[];
for i=1:0.1:4
    % n_f2_calculados=linspace(2.415*i,4*2.415*i,4);
    % for j=n_f2_calculados
    for j=10:-0.1:2.415*i
        freq_list = [freq_list ; [i,j]];
    end
end


%%
zeta1=0.1;
zeta2=0.25;

for i = 1:size(freq_list, 1)
    f1 = freq_list( i, 1);
    f2 = freq_list( i, 2);

    % --- Controller and transfer function computations ---
    k_p=1.2993/1e-2;
    G_c = tf(k_p,1);
    [~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,~,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , mass,mass , f1, zeta1 , f2 , zeta2);
    tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
    G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
    [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , mass,mass , f1, zeta1 , f2 , zeta2);

    x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
    ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');

    % i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh'); % converting mm to m
    % F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');

    [picos_ddx_table , picos_x_table ] = ResponseSpectrum( t_vector , ddx_T, f_vector );

    % Compute MSE for acceleration response (you can also compute for displacement)
    mse = mean((picos_ddx_table - picos_ddx_ground).^2);

    % Save the values for the 3D plot (here we use the acceleration MSE)
    f1_arr(end+1) = f1;
    f2_arr(end+1) = f2;
    mse_arr(end+1)   = mse;
    % isv_arr(end+1) = max(abs(i_sv));
    % Fp_arr(end+1)= max(abs(F_p_isv));
end


%% Generate a 3D plot: MSE vs. f1 and f2
figure(1);
f1_unique = unique(f1_arr);
f2_unique = unique(f2_arr);
[Z1, Z2] = meshgrid(f1_unique, f2_unique);

% Interpolating scattered data onto a grid
MSE_matrix = griddata(f1_arr, f2_arr, mse_arr, Z1, Z2, 'linear');

% Surface plot
surf(Z1, Z2, MSE_matrix);
xlabel('f_1');
ylabel('f_2');
zlabel('MSE');
title(sprintf('3D Plot of MSE vs. f_1 and f_2 (m_i=%.0f ton ,ξ_1=%.0f Hz & ξ_2=%.0f Hz)', mass*1e-3,zeta1,zeta2));
grid on;
colormap(jet(256));
colorbar;
saveas(fig1, sprintf('3D Plot of MSE vs. f_1 and f_2 (m_i=%.0f ton ,ξ_1=%.0f Hz & ξ_2=%.0f Hz)', mass*1e-3,zeta1,zeta2);
 
