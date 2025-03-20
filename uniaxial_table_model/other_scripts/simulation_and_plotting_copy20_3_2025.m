clear 
close all
%clc

%% Structure parameters

mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;

% 1st mode
m1 = mass; % kg
f1 = 0.4; % Hz   % 1.5 < f1 < 4
zeta1 = 0.1 ; % 2 < zeta1 < 10
%2nd mode
m2 = mass; % kg
f2 =3; % Hz % 6 < f2 < 10
zeta2 = 0.06; % 5 < zeta2 < 25

% Controller
k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2, ss_model 
[s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ss_isv_xT  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2); %_and_StateSpace
%ss_xref_xT = feedback(ss_isv_xT*G_c ,1);

%%  Load seismic signal and scale down if necessary
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2);
ddx_ref = dados(:,2);
ddx = [t_vector dados(:,2)];
%ddy = [t_vector  dados(:,3)];

% Limits
lim_displacement = 100; % mm
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

% Scaling down if necessary
% displacements in milimeters
x_ref = lsim(1e3/s^2,  ddx_ref , t_vector ,'foh');
max_xref = max(x_ref);
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

%% Plots

fig1 = figure(1); %clf; % Clear figure if needed
ax1 = axes(fig1); hold(ax1, 'on');
opts1=bodeoptions('cstprefs');
opts1.FreqUnits = 'Hz';
opts1.XLim={[1 40]};
% opts1.YLim={[-40 1]};
% opts1.MagVisible='off';

fig2 = figure(2); %clf;
ax2 = axes(fig2); hold(ax2, 'on');
grid on
hold on
title('Input to Servo'); 
legend();
xlabel('Time (s)'); 
ylabel('Voltage (V)');

fig3 = figure(3); %clf;
ax3 = axes(fig3); hold(ax3, 'on');
grid on
title(' Platen Displacement '); 
hold on
legend()
xlabel('Time (s)'); 
ylabel('Displacement (mm)');

fig4 = figure(4); %clf;
ax4 = axes(fig4); hold(ax4, 'on');
grid on
title('Platen Acceleration'); 
hold on
legend()
xlabel('Time (s)'); 
ylabel('Acceleration (m/s^2)');

fig5= figure(5); %clf; 
ax5 = axes(fig5); hold(ax5, 'on');
grid on
title(' Platen Displacement Tracking Error '); 
legend()
xlabel('Time (s)'); 
ylabel('Error (mm)');

fig6= figure(6); %clf; 
ax6 = axes(fig6); hold(ax6, 'on');
grid on
title(' Platen Acceleration Tracking Error '); 
legend()
xlabel('Time (s)'); 
ylabel('Error (m/s^2)');

fig7= figure(7); %clf; 
ax7 = axes(fig7); hold(ax7, 'on');
grid on
title('Force to Platen'); 
legend();
xlabel('Time (s)'); 
ylabel('Force (kN)');

fig8 = figure(8);%
subplot(121);
grid on;
xlabel('Frequency (Hz)');
ylabel('Acceleration (m/s^2)');
title('Acceleration Response Spectra');
xlim([1 30]);
subplot(122);
grid on;
xlabel('Frequency (Hz)');
ylabel('Displacement (m)');
title('Displacement Response Spectra');
xlim([0.1 5]);
% Define colors for lines 1/3 and 2/4
color1 = 'r'; % MATLAB default blue
color2 = 'b'; % MATLAB default orange
color3 = 'y';


%% Finding Response Spectre of Ground

f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , ddx_ref, f_vector , 1);

figure(fig8);
subplot(121)
grid on;
legend();
hold on
plot(f_vector, picos_ddx_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground');% - Normal
%semilogx(f_vector, picos_ddx_ground(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Ground - Parallel');

subplot(122)
grid on;
legend();
hold on
plot(f_vector, picos_x_ground(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Ground ');%- Normal
%semilogx(f_vector, picos_x_ground(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Ground - Parallel');


%%
% First plot
axes(ax1); % Activate the existing axes
bodeplot(G_xT_xref,opts1);
% hold on
% bodeplot(ss_xref_xT)

% Third plot
axes(ax3); % Activate the existing axes

x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
erro = x_T-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_ref,"DisplayName","Reference")
plot(t_vector,x_T,"DisplayName","MSE="+string(mse))
%x_T_ss = lsim(ss_xref_xT,x_ref,t_vector);
%plot(t_vector,x_T_ss,'--',"DisplayName","State Space")

%%
% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');
erro = ddx_T-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_ref,"DisplayName","Reference")
plot(t_vector,ddx_T,"DisplayName","MSE="+string(mse))

% 6th plot
axes(ax6); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Default")

% Second plot
axes(ax2); % Activate the existing axes
i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","Default")

% 7th  plot 
axes(ax7); % Activate the existing axes
F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
%plot(t_vector,ddx_ref*mT/1e3,"DisplayName","Reference") % weird values
plot(t_vector,F_p_isv/1e3,"DisplayName","Default") %/1e3 to display as kN



%% Finding Response Spectre for table

[picos_ddx_table , picos_x_table ] = ResponseSpectrum(  t_vector , ddx_T , f_vector, 1 );

figure(fig8);
subplot(121)
hold on
mse = mean((picos_ddx_table-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName',  sprintf('Platform - MSE= %.2e', mse(1)));
%semilogx(f_vector, picos_ddx_table(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', "MSE="+string(mse(2)));%- Parallel

subplot(122)
hold on
mse = mean((picos_x_table-picos_x_ground).^2);
plot(f_vector, picos_x_table(:, 1),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', sprintf('Platform - MSE= %.2e', mse(1)));
%semilogx(f_vector, picos_x_table(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', "MSE="+string(mse(2)));



%% Use PID tuner app to generate a PID controller for the system
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
[s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ ,~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

%%
% First plot
axes(ax1); % Activate the existing axes
bodeplot(G_xT_xref);
legend( 'Default'  ,  'Tuned');
title('Bode of G\_xT\_xref'); 
grid on;

% displacements in milimeters
axes(ax3); % Activate the existing axes
x_T_tuned = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
erro = x_T_tuned-x_ref;
mse = mean(erro.^2);
plot(t_vector,x_T_tuned,"DisplayName","Tuned MSE="+string(mse))

% 5th plot
axes(ax5); % Activate the existing axes
plot(t_vector,erro,"DisplayName","Tuned")

% Fourth plot
axes(ax4); % Activate the existing axes
ddx_T_tuned = lsim(G_xT_xref, ddx_ref ,t_vector,'foh');
erro = ddx_T_tuned-ddx_ref;
mse = mean(erro.^2);
plot(t_vector,ddx_T_tuned,"DisplayName","Tuned MSE="+string(mse))

% 6th plot
axes(ax6); % Activate the existing axes
hold on
plot(t_vector,erro,"DisplayName","Tuned")


% Second plot
axes(ax2); % Activate the existing axes
hold on
i_sv = lsim(G_c ,   (x_ref-x_T_tuned)*1e-3  ,t_vector,'foh');
plot(t_vector,i_sv,"DisplayName","Tuned")
%  plot
axes(ax7); % Activate the existing axes
hold on
F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
plot(t_vector,F_p_isv*1e-3,"DisplayName","Tuned")



%% Finding Response Spectre for table tuned

[picos_ddx_table_tuned , picos_x_table_tuned ] = ResponseSpectrum( t_vector , ddx_T_tuned, f_vector , 1 );


figure(fig8);
subplot(121)
hold on
mse = mean((picos_ddx_table_tuned-picos_ddx_ground).^2);
plot(f_vector, picos_ddx_table_tuned(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName', sprintf('Tuned Platform - MSE= %.2e', mse(1)));
%semilogx(f_vector, picos_ddx_table_tuned(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Tuned Platform - Parallel');

subplot(122)
hold on
mse = mean((picos_x_table_tuned-picos_x_ground).^2);
plot(f_vector, picos_x_table_tuned(:, 1),'-', 'LineWidth' , 2, 'Color', color3, 'DisplayName',  sprintf('Tuned Platform - MSE= %.2e', mse(1)));
%semilogx(f_vector, picos_x_table_tuned(:, 2),'-', 'LineWidth' , 2, 'Color', color2, 'DisplayName', 'Tuned Platform - Parallel');




%% State Space model




% % u = Fp
% % y = x_ss
% 
% Css= eye(6);
% Dss = 0;
% 
% MV = struct(Min=-lim_force,Max=lim_force);
% p = 20;
% m = 3;
% mpcobj = mpc(Plant,Ts,p,m,[],MV);
% 



%% Save all figures after plotting
% Folder path where you want to save the images
folderName = sprintf('Sim_Res/m_i=%.1f,f_1=%.1f, f_2=%.1f,zeta1=%.2f,zeta2=%.2f',mass*1e-3,f1,f2,zeta1,zeta2);

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % Create the folder if it doesn't exist
    mkdir(folderName);
end

saveas(fig1,fullfile(folderName,'Bode_of_G_xT_xref.png'));
saveas(fig2,fullfile(folderName,'Input_to_Servo.png'));
saveas(fig3,fullfile(folderName,'Platen_Displacement.png'));
saveas(fig4,fullfile(folderName,'Platen_Acceleration.png'));
saveas(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.png'));
saveas(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.png'));
saveas(fig7,fullfile(folderName,'Force_to_Platen.png'));
saveas(fig8,fullfile(folderName,'Response_Spectra.png'));
 
saveas(fig1,fullfile(folderName,'Bode_of_G_xT_xref.fig'));
saveas(fig2,fullfile(folderName,'Input_to_Servo.fig'));
saveas(fig3,fullfile(folderName,'Platen_Displacement.fig'));
saveas(fig4,fullfile(folderName,'Platen_Acceleration.fig'));
saveas(fig5,fullfile(folderName,'Platen_Displacement_Tracking_Error.fig'));
saveas(fig6,fullfile(folderName,'Platen_Acceleration_Tracking_Error.fig'));
saveas(fig7,fullfile(folderName,'Force_to_Platen.fig'));
saveas(fig8,fullfile(folderName,'Response_Spectra.fig'));
