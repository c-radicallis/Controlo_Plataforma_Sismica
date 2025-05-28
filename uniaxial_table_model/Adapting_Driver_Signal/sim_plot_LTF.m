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
[~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , AA , BB , CC , DD  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

%% --- Simulation --

addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'\Adapting_Driver_Signal\LTF_to_TXT\

loadLNEC('LAquilaReducedScale.tgt.txt')

loadLNEC('LAquilaReducedScale_34.DRV.txt')

filename = 'TestSequence.xlsx'; % Define the Excel file and sheet
sheet = 1;  % or use sheet name, e.g., 'Sheet1'
data = readtable(filename, 'Sheet', sheet);% Read the Excel file into a table
targetDRV = 'LAquilaReducedScale_34.DRV';  % Define the target DRV value to search for
rowIndex = find(strcmp(data.DRV, targetDRV), 1);% Find the row index where DRV matches the target
if isempty(rowIndex)
    disp('DRV value not found.');
else    
    scaleFactor = data.ScaleFactor(rowIndex)
end

plot(time_vector , scaleFactor*tgtDispT)
hold on
plot(time_vector , DRVDispT )

%%
opts = detectImportOptions('LTF_to_TXT/LAquilaReducedScale_34_DRV.txt', 'FileType','text');
opts.DataLines = [2 Inf];  % tell it that the first line is header:
dados = readmatrix('LTF_to_TXT/LAquilaReducedScale_34_DRV.txt', opts);
t_vector = dados(:,1);
t_step = t_vector(2);
x_drv    = dados(:,2);



opts = detectImportOptions('LTF_to_TXT\LAquilaReducedScale_tgt.txt', 'FileType','text');
opts.DataLines = [2 Inf];  % tell it that the first line is header:
dados = readmatrix('LTF_to_TXT\LAquilaReducedScale_tgt.txt', opts);
dados = scaleFactor*dados;
x_tgt = dados(:,2);
ddx_tgt = dados(:,3);

%%
plot(t_vector,x_tgt);hold on;
double_int_x_tgt = cumtrapz(cumtrapz(x_tgt,t_vector),t_vector);
plot(t_vector,double_int_x_tgt)

%% Plots
close all;
fig1 = figure(1); ax1 = axes(fig1); hold(ax1, 'on'); opts1=bodeoptions('cstprefs'); opts1.FreqUnits = 'Hz'; opts1.XLim={[1 50]}; % opts1.YLim={[-40 1]}; % opts1.MagVisible='off';
fig2 = figure(2); ax2 = axes(fig2); hold(ax2, 'on'); grid on; hold on; title('Input to Servo');  legend(); xlabel('Time (s)');  ylabel('Voltage (V)');
fig3 = figure(3); ax3 = axes(fig3); hold(ax3, 'on'); grid on; title(' Platen Displacement ');  hold on; legend(); xlabel('Time (s)');  ylabel('Displacement (mm)');
fig4 = figure(4); ax4 = axes(fig4); hold(ax4, 'on'); grid on; title('Platen Acceleration');  hold on; legend(); xlabel('Time (s)');  ylabel('Acceleration (m/s^2)'); 
fig5 = figure(5); ax5 = axes(fig5); hold(ax5, 'on');grid on;title(' Platen Displacement Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (mm)');
fig6 = figure(6); ax6 = axes(fig6); hold(ax6, 'on');grid on;title(' Platen Acceleration Tracking Error '); legend();xlabel('Time (s)'); ylabel('Error (m/s^2)');
fig7 = figure(7); ax7 = axes(fig7); hold(ax7, 'on');grid on;title('Force to Platen'); legend();xlabel('Time (s)'); ylabel('Force (kN)');
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra');xlim([1 30]);subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; % Define colors for lines 1/3 and 2/4

%% Finding Response Spectre of Ground
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_tgt , picos_x_tgt] = ResponseSpectrum( t_vector , ddx_tgt, f_vector , 1);

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal

%% Use PID tuner app to generate a PID controller for the system
axes(ax1);
bodeplot(G_xT_xref);

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
[s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref_tuned,~,~ , G_Fp_isv  ,~,~,~,~ ,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);
x_T_tuned = lsim(G_xT_xref_tuned/s^2 ,  ddx_tgt ,t_vector,'foh');
ddx_T_tuned = lsim(G_xT_xref_tuned, ddx_tgt ,t_vector,'foh');

axes(ax1); hold on;
bodeplot(G_xT_xref_tuned);

axes(ax3); % Activate the existing axes
plot(t_vector,x_tgt,"DisplayName","Reference")
erro = x_T_tuned-x_tgt;
mse = mean(erro.^2);
plot(t_vector,x_T_tuned,"DisplayName","Tuned MSE="+string(mse))

axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Tuned")

axes(ax4); % Activate the existing 
plot(t_vector,ddx_tgt,"DisplayName","Reference")
erro = ddx_T_tuned-ddx_tgt;
mse = mean(erro.^2);
plot(t_vector,ddx_T_tuned,"DisplayName","Tuned MSE="+string(mse))

axes(ax6); hold on;
plot(t_vector,erro,"DisplayName","Tuned")

axes(ax2); hold on;
i_sv_tuned = lsim(G_c ,  x_tgt-x_T_tuned  ,t_vector,'foh');
plot(t_vector,i_sv_tuned,"DisplayName","Tuned")

axes(ax7); hold on;
F_p_isv = lsim(G_Fp_isv,   i_sv_tuned  , t_vector,'foh');
plot(t_vector,F_p_isv*1e-3,"DisplayName","Tuned")

%% Finding Response Spectre for table tuned
[picos_ddx_table_tuned , picos_x_table_tuned ] = ResponseSpectrum( t_vector , ddx_T_tuned, f_vector , 1 );

figure(fig8);subplot(121);hold on;
mse = mean((picos_ddx_table_tuned-picos_ddx_tgt).^2);
plot(f_vector, picos_ddx_table_tuned(:, 1),'-', 'LineWidth' , 2,  'DisplayName', sprintf('Tuned Platform - MSE= %.2e', mse(1)));
subplot(122);hold on;
mse = mean((picos_x_table_tuned-picos_x_tgt).^2);
plot(f_vector, picos_x_table_tuned(:, 1),'-', 'LineWidth' , 2,  'DisplayName',  sprintf('Tuned Platform - MSE= %.2e', mse(1)));

%% Optimal control
sys = ss(AA,BB,CC,DD);
nx = size(AA,1);    % Number of states
nu = size(BB,2);    % Number of control inputs (should be 1)
ny = size(CC,1);    % Number of outputs

plant_aug = ss(AA, BB,[eye(nx);CC],DD);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'y_xT'};  % plant output

sumblk1 = sumblk('e = x_tgt - y_xT'); % Compute the error signal: e = r - y

integrator = tf(1,[1 0]); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = 1e3*diag([zeros(1,nx),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-9*eye(nu);
K_lqi = lqi(sys, Q, R)% Design the LQI controller for the original system
K  = K_lqi(1:nx);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = {'Qsv' , 'Fp' , 'xT' , 'x1' , 'x2','dxT' , 'dx1' , 'dx2' , 'xi'};
controller.OutputName = {'i_sv'};

clsys = connect(plant_aug,  controller , integrator, sumblk1, 'x_tgt', 'y_xT')

[x_T_LQG, t_out, x] = lsim(clsys, x_tgt, t_vector,'foh'); % Simulate the closed-loop response using lsim:
diff_x_T_LQG=diff(x_T_LQG)./diff(t_out);
ddx_T_LQG=diff(diff_x_T_LQG)./diff(t_out(1:end-1));
%% Ploting LQG results
axes(ax1);
bodeplot(clsys,opts1);
legend( 'Default' , 'PIDF'  ,'LQG');
title('Bode of G\_xT\_xref'); 
grid on;

axes(ax3);% Third plot % Activate the existing axes
erro = x_T_LQG-x_tgt;
mse = mean(erro.^2);
plot(t_vector,x_T_LQG,"DisplayName","MSE="+string(mse))

axes(ax5); hold on;% 5th plot % Activate the existing axes
plot(t_vector,erro,"DisplayName","LQG")

axes(ax4); hold on; % Activate the existing axes
erro = ddx_T_LQG-ddx_tgt(1:end-2);
mse = mean(erro.^2);
plot(t_vector(1:end-2),ddx_T_LQG,"DisplayName","MSE="+string(mse))

axes(ax6); hold on;% Activate the existing axes
plot(t_vector(1:end-2),erro,"DisplayName","LQG")

axes(ax2); hold on;
i_sv = lsim(clsys ,   x_tgt - x_T_LQG  ,t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","LQG")

axes(ax7); hold on;
F_p_isv = x(:,2); % 2nd element of state vector
plot(t_vector,F_p_isv*1e-3,"DisplayName","LQG")

%% Finding Response Spectre for table LQG
[picos_ddx_table_LQG , picos_x_table_LQG ] = ResponseSpectrum( t_vector , ddx_T_LQG, f_vector , 1 );

figure(fig8); subplot(121); hold on;
mse = mean((picos_ddx_table_LQG-picos_ddx_tgt).^2);
plot(f_vector, picos_ddx_table_LQG(:, 1),'-', 'LineWidth' , 2,  'DisplayName', sprintf('LQG Platform - MSE= %.2e', mse(1)));

subplot(122);hold on;
mse = mean((picos_x_table_LQG-picos_x_tgt).^2);
plot(f_vector, picos_x_table_LQG(:, 1),'-', 'LineWidth' , 2,  'DisplayName',  sprintf('LQG Platform - MSE= %.2e', mse(1)));

%% Save all figures after plotting
folderName = sprintf('Sim_Res_LQG/Q=%.1e,m_i=%.1f,f_1=%.1f, f_2=%.1f,zeta1=%.2f,zeta2=%.2f',Q(9,9),mass*1e-3,f1,f2,zeta1,zeta2); % Folder path where you want to save the images
if ~exist(folderName, 'dir')% Check if the folder already exists
    % Create the folder if it doesn't exist
    mkdir(folderName);
end

drawnow expose

pause(10)

tightfig(fig2); tightfig(fig3); tightfig(fig4); tightfig(fig5); tightfig(fig6); tightfig(fig7); tightfig(fig8); 


exportgraphics(fig1,fullfile(folderName,'Bode_of_G_xT_xref.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig2,fullfile(folderName,'Input_to_Servo.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig3,fullfile(folderName,'Platen_Displacement.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig4,fullfile(folderName,'Platen_Acceleration.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig7,fullfile(folderName,'Force_to_Platen.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
exportgraphics(fig8,fullfile(folderName,'Response_Spectra.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
 
saveas(fig1,fullfile(folderName,'Bode_of_G_xT_xref.fig'));
saveas(fig2,fullfile(folderName,'Input_to_Servo.fig'));
saveas(fig3,fullfile(folderName,'Platen_Displacement.fig'));
saveas(fig4,fullfile(folderName,'Platen_Acceleration.fig'));
saveas(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.fig'));
saveas(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.fig'));
saveas(fig7,fullfile(folderName,'Force_to_Platen.fig'));
saveas(fig8,fullfile(folderName,'Response_Spectra.fig'));

