clear 
close all
%clc

%%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);

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
f_i=0.01; %freq inicial
f_n=30;  %freq final
n_points = 200;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

[picos_ddx_ground , ~] = ResponseSpectrum( t_vector , ddx_ref , f_vector , 0);

%% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;


%%
zeta_list =[];
% for i=1:0.2:6
%     for j=4:24
% for i=1:0.1:6
%     for j=4:0.5:24
for i=1:15
    for j=1:25
        zeta_list = [zeta_list ; [i,j]];
    end
end
zeta_list = zeta_list/100;

%%

f1 = 0.4;
f2 = 3.5;

elements =size(zeta_list, 1);
% Before the loops, initialize arrays to store zeta1, zeta2, and mse values
zeta1_arr = zeros(elements,1);
zeta2_arr =zeros(elements,1);
mse_arr   = zeros(elements,1);
% isv_arr=zeros(elements,1);
% Fp_arr=zeros(elements,1);

for i= 1:elements
    sprintf('Step %.0f of %.0f',i,elements)
    zeta1 = zeta_list( i, 1);
    zeta2 = zeta_list( i, 2);

    % --- Controller and transfer function computations ---
    k_p=1.2993/1e-2;  %V/m
    G_c = tf(k_p,1);
    [~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,~,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , mass,mass , f1, zeta1 , f2 , zeta2);
    tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
    G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
    [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , mass,mass , f1, zeta1 , f2 , zeta2);

    %x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
    ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');

    % i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh'); % converting mm to m
    % F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
    % isv_arr(end+1) = max(abs(i_sv));
    % Fp_arr(end+1)= max(abs(F_p_isv));

    [picos_ddx_table , ~ ] = ResponseSpectrum( t_vector , ddx_T, f_vector , 0);

    % Compute MSE for acceleration response (you can also compute for displacement)
    mse = mean((picos_ddx_table - picos_ddx_ground).^2);

    % Save the values for the 3D plot (here we use the acceleration MSE)
    zeta1_arr(i) = zeta1;
    zeta2_arr(i) = zeta2;
    mse_arr(i)   = mse;

end


%% Scatter Generate a 3D plot: MSE vs. zeta1 and zeta2
fig1 = figure(1);
[Z1, Z2] = meshgrid(unique(zeta1_arr), unique(zeta2_arr));
MSE_matrix = reshape(mse_arr, size(Z1));
surf(Z1, Z2, MSE_matrix)%,"FaceColor","interp");%,"EdgeColor","interp"
%scatter3(zeta1_arr, zeta2_arr, mse_arr, 100, mse_arr, 'filled');
xlabel('\zeta_1');
ylabel('\zeta_2');
zlabel('MSE');
title(sprintf('3D Plot of MSE vs. ξ_1 and ξ_2 (m_i=%.1f ton , f_1=%.1f Hz & f_2=%.1f Hz)', mass*1e-3,f1,f2));
grid on;
colormap(jet(256));
colorbar;

% Folder path where you want to save the images
folderName = 'MSE_Zeta';

% % Check if the folder already exists
% if ~exist(folderName, 'dir')
%     % Create the folder if it doesn't exist
%     mkdir(folderName);
% end

%%

saveas(fig1, fullfile(folderName,sprintf('MSE vs Zeta (m_i=%.1f,f_1=%.1f, f_2=%.1f).fig', mass*1e-3,f1,f2)));
%exportgraphics(gca,fullfile(folderName,sprintf('MSE vs Zeta (m_i=%.1f,f_1=%.1f, f_2=%.1f).fig', mass*1e-3,f1,f2)),'BackgroundColor','none');
