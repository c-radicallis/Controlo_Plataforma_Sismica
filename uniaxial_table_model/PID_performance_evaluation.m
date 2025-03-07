clear 
close all
%clc

%% Plots

fig2 = figure(2); %clf;
ax2 = axes(fig2); hold(ax2, 'on');
grid on
hold on
title('Input to Servo'); 
legend();
xlabel('Time (s)'); 
ylabel('Voltage (V)');

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
ylabel('ddx_m (m/s^2)');
title('Response Spectra');
subplot(122);
grid on;
xlabel('Frequency (Hz)');
ylabel('x_m (m)');
title('Response Spectra');
xlim([0 5]);
% Define colors for lines 1/3 and 2/4
color1 = 'r'; % MATLAB default blue
color2 = 'b'; % MATLAB default orange
color3 = 'g';

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

%% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=[1 , 2]*1e3;
for m_i = mass

    % 1st mode
    m1 = mass; % kg
    f1 = 2; % Hz   % 1.5 < f1 < 4
    zeta1 = 0.05 ; % 2 < zeta1 < 10
    %2nd mode
    m2 = m1; % kg
    f2 = 10; % Hz % 6 < f2 < 10
    zeta2 = 0.1; % 5 < zeta2 < 25

    % Use PID tuner app to generate a PID controller for the system
    tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
    G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
    [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

    % Finding Response Spectre of Ground

    f_i=0.1; %freq inicial
    f_n=30;  %freq final
    n_points = 1e2;
    f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
    [filtered_picos_ddx_ground , filtered_picos_x_ground] = ResponseSpectre_filtered( ddx_ref, f_vector );

    figure(fig8);
    subplot(121)
    grid on;
    legend();
    hold on
    semilogx(f_vector, filtered_picos_ddx_ground(:, 1),'-o', 'LineWidth' , 1, 'Color', color1, 'DisplayName', 'Ground');% - Normal

    subplot(122)
    grid on;
    legend();
    hold on
    semilogx(f_vector, filtered_picos_x_ground(:, 1),'-o', 'LineWidth' , 1, 'Color', color1, 'DisplayName', 'Ground ');%- Normal


    %

    x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
    ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');

    % Second plot
    axes(ax2); % Activate the existing axes
    i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh'); % 1e-3 converts from mm to meters
    plot(t_vector,i_sv,"DisplayName","Default")

    % 7th  plot 
    axes(ax7); % Activate the existing axes
    F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
    plot(t_vector,F_p_isv*1e-3,"DisplayName","Default") % 1e-3 to display as kN


    % Finding Response Spectre for table
    [filtered_picos_ddx_table , filtered_picos_x_table ] = ResponseSpectre_filtered( [ t_vector , ddx_T], f_vector );

    figure(fig8);
    subplot(121)
    hold on
    mse = mean((filtered_picos_ddx_table-filtered_picos_ddx_ground).^2);
    semilogx(f_vector, filtered_picos_ddx_table(:, 1),'-+', 'LineWidth' , 1, 'Color', color2, 'DisplayName', "Platform - MSE="+string(mse(1))); % - Normal

    subplot(122)
    hold on
    mse = mean((filtered_picos_x_table-filtered_picos_x_ground).^2);
    semilogx(f_vector, filtered_picos_x_table(:, 1),'-+', 'LineWidth' , 1, 'Color', color2, 'DisplayName',"Platform - MSE="+string(mse(1)));


end

%% Save all figures after plotting
saveas(fig2, 'Input_to_Servo.png');
saveas(fig7, 'Force_to_Platen.png');
saveas(fig8, 'Response_Spectra.png');


%%

function [ filtered_picos_ddx_m , filtered_picos_x_m    ] = ResponseSpectre_filtered( dados , f_vector )

     t_vector=dados(:,1);
    [ ~ , n_directions]=size(dados);

    m=1;%1kg
    zeta=0.05; %damping ratio (%)
    
    picos_ddx_m = zeros( length(f_vector)  , n_directions-1 );
    picos_x_m = picos_ddx_m;
    
    for i=1:length(f_vector)    
      
        % from the desired natural frequency, determining stiffness and damping  
        k = m*(2*pi*f_vector(i))^2; %N/m
        c = zeta*2*m*2*pi*f_vector(i); %N/m/s
          
        for j=2:n_directions
            % simulacao do modelo 
            s=tf('s');
            ddx_m = lsim( (c*s+k)/(m*s^2+c*s+k) , dados(:,j) , t_vector ,'zoh'); % Funçao de tranferencia de mola massa amortecedor
            picos_ddx_m(i,j-1)=max(abs( ddx_m(:,1) ));  
    
            x_m = lsim( (c*s+k)/(m*s^2+c*s+k)*1/s^2 , dados(:,j)  , t_vector ,'zoh');
            picos_x_m(i,j-1)=max(abs( x_m(:,1) ));  
    
        end
    end

    % Apply a filter to the data
    % % Moving average filter
    windowSize = length(f_vector)*0.3; % Adjust window size as needed as percentage of elements in f_vector
    filtered_picos_ddx_m             = movmean(picos_ddx_m, windowSize);
    filtered_picos_x_m                 = movmean(picos_x_m, windowSize);

end
