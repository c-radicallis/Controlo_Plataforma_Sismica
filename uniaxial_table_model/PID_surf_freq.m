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

[filtered_picos_ddx_ground , filtered_picos_x_ground] = ResponseSpectre_filtered( t_vector , ddx_ref , f_vector );

%% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)
m1=2e3;
m2=m1;

% Before the loops, initialize arrays to store f1, f2, and mse values
f1_arr = [];
f2_arr = [];
mse_arr   = [];
isv_arr=[];
Fp_arr=[];

%%
freq_list =[];
for i=1:4
    for j=2.415*i:2.415:4*2.415*i
        freq_list = [freq_list ; [i,j]];
    end
end

%%
zeta1=0.02;
zeta2=0.05;

for i = 1:size(freq_list, 1)
    f1 = freq_list( i, 1);
    f2 = freq_list( i, 2);

    % --- Controller and transfer function computations ---
    k_p=1.2993/1e-2;
    G_c = tf(k_p,1);
    [~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,~,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);
    tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
    G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
    [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

    x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
    ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');

    % i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh'); % converting mm to m
    % F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');

    [filtered_picos_ddx_table , filtered_picos_x_table ] = ResponseSpectre_filtered( t_vector , ddx_T, f_vector );

    % Compute MSE for acceleration response (you can also compute for displacement)
    mse = mean((filtered_picos_ddx_table - filtered_picos_ddx_ground).^2);

    % Save the values for the 3D plot (here we use the acceleration MSE)
    f1_arr(end+1) = f1;
    f2_arr(end+1) = f2;
    mse_arr(end+1)   = mse;
    % isv_arr(end+1) = max(abs(i_sv));
    % Fp_arr(end+1)= max(abs(F_p_isv));
end


%% Generate a 3D plot: MSE vs. f1 and f2
figure(1);
[Z1, Z2] = meshgrid(unique(f1_arr), unique(f2_arr));
MSE_matrix = reshape(mse_arr, size(Z1));
surf(Z1, Z2, MSE_matrix);

xlabel('f_1');
ylabel('f_2');
zlabel('MSE');
title('3D Plot of MSE vs. f_1 and f_2');
grid on;
colorbar;

%% Funtions
function [ filtered_picos_ddx_m , filtered_picos_x_m    ] = ResponseSpectre_filtered( t_vector , accel , f_vector )

    m=1;%1kg
    zeta=0.05; %damping ratio (%)s
    
    picos_ddx_m = zeros( length(f_vector) ,1 );
    picos_x_m = picos_ddx_m;
    
    for i=1:length(f_vector)    
      
        % from the desired natural frequency, determining stiffness and damping  
        k = m*(2*pi*f_vector(i))^2; %N/m
        c = zeta*2*m*2*pi*f_vector(i); %N/m/s
          
        % simulacao do modelo 
        s=tf('s');
        ddx_m = lsim( (c*s+k)/(m*s^2+c*s+k) , accel , t_vector ,'zoh'); % Funçao de tranferencia de mola massa amortecedor
        picos_ddx_m(i)=max(abs( ddx_m(:,1) ));  

        x_m = lsim( (c*s+k)/(m*s^2+c*s+k)*1/s^2 , accel  , t_vector ,'zoh');
        picos_x_m(i)=max(abs( x_m(:,1) ));  

    end

    % Apply a filter to the data
    % % Moving average filter
    windowSize = length(f_vector)*0.15; % Adjust window size as needed as percentage of elements in f_vector
    filtered_picos_ddx_m = movmean(picos_ddx_m, windowSize);
    filtered_picos_x_m  = movmean(picos_x_m, windowSize);

end

