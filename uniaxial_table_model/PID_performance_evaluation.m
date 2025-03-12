clear 
close all
%clc

%% Plots

% fig2 = figure(2); %clf;
% ax2 = axes(fig2); hold(ax2, 'on');
% grid on
% hold on
% title('Input to Servo'); 
% legend();
% xlabel('Time (s)'); 
% ylabel('Voltage (V)');
% 
% fig7= figure(7); %clf; 
% ax7 = axes(fig7); hold(ax7, 'on');
% grid on
% title('Force to Platen'); 
% legend();
% xlabel('Time (s)'); 
% ylabel('Force (kN)');

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
% % Define colors for lines 1/3 and 2/4
% color1 = 'r'; % MATLAB default blue
% color2 = 'b'; % MATLAB default orange
% color3 = 'g';

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
n_points = 500;
%f_vector = linspace( f_i , f_n , n_points); 
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

[filtered_picos_ddx_ground , filtered_picos_x_ground] = ResponseSpectre_filtered( t_vector , ddx_ref , f_vector );

%% Structure parameters
mT=1.9751*1e3;       %Platen mass (mp=1.9751 t)
cT=5.78*1e3;        %Total damping, actuator + platen (ct=5.78 kN s/m1)

mass=2e3;
for m_i = mass
    m1 = m_i; % kg
    m2 = m_i; % kg

    % 1.5 < f1 < 4
    % 6 < f2 < 10
    f_list =  [ [1.5,6] ];% ; [4,10] ]; % Define the list
    for i = 1:size(f_list, 1)
        f1 = f_list( i, 1);
        f2 = f_list( i, 2);
        
        % 2 < zeta1 < 10
        % 5 < zeta2 < 25
        zeta_list = [ [2 , 5] ; [2 , 25] ; [10 , 25] ; [2 , 25]]/100; % Define the list
        for j = 1:size(zeta_list, 1)
            zeta1 = zeta_list( j, 1);
            zeta2 = zeta_list( j, 2);

            % Controller
            k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
            G_c = tf(k_p,1);
            [~,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,~,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);
            % Use PID tuner app to generate a PID controller for the system
            tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
            G_c   = pidtune(G_Fp_isv*G_xT_Fp,'PIDF',20*2*pi,tuner_opts)
            [s,~,~,~,~ ,~ ,~,~,~,~,G_xT_Fp,~,G_xT_xref,~,~ , G_Fp_isv  ,~,~,~,~ , ~ ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

            %

            x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx_ref ,t_vector,'foh');
            ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');

            % % Second plot
            % axes(ax2); % Activate the existing axes
            i_sv = lsim(G_c,   (x_ref-x_T)*1e-3  , t_vector,'foh'); % 1e-3 converts from mm to meters
            max_isv = max(abs(i_sv));
            % plot(t_vector,i_sv,"DisplayName", sprintf('m=%i , f_1=%.1f , f_2=%.1f  , ξ_1=%.2f ,  ξ_2=%.2f',m_i , f1, f2 , zeta1 , zeta2 ))
            % 
            % % 7th  plot 
            % axes(ax7); % Activate the existing axes
            F_p_isv = lsim(G_Fp_isv,   i_sv  , t_vector,'foh');
            max_Fp = max(abs(F_p_isv));
            % plot(t_vector,F_p_isv*1e-3,"DisplayName",sprintf('m=%i , f_1=%.1f , f_2=%.1f  , ξ_1=%.2f ,  ξ_2=%.2f',m_i , f1, f2 , zeta1 , zeta2 )) % 1e-3 to display as kN


            % Finding Response Spectre for table
            [filtered_picos_ddx_table , filtered_picos_x_table ] = ResponseSpectre_filtered( t_vector , ddx_T, f_vector );

            figure(fig8);
            subplot(121)
            hold on
            mse = mean((filtered_picos_ddx_table-filtered_picos_ddx_ground).^2);
            semilogx(f_vector, filtered_picos_ddx_table,'-', 'LineWidth' , 1,  'DisplayName', sprintf(' f_{1}=%.1f , f_{2}=%.1f , ξ_{1}=%.2f  , ξ_{2}=%.2f  , i_{sv}= %.1f , F_{P}=%.0f , MSE=%.1e' , f1, f2 , zeta1 , zeta2 ,max_isv , round(max_Fp*1e-3,3) , mse) ); % - Normal

            subplot(122)
            hold on
            mse = mean((filtered_picos_x_table-filtered_picos_x_ground).^2);
            semilogx(f_vector, filtered_picos_x_table,'-', 'LineWidth' , 1,  'DisplayName', sprintf(' f_{1}=%.1f , f_{2}=%.1f  , ξ_{1}=%.2f ,  ξ_{2}=%.2f , i_{sv}= %.1f , F_{P}=%.0f , MSE=%.1e' , f1, f2 , zeta1 , zeta2 ,max_isv , round(max_Fp*1e-3,3) , mse));
            % sprintf(' f_1=%.1f , f_2=%.1f  , ξ_1=%.2f ,  ξ_2=%.2f  i_{sv}= %.1f  , F_P=%.1e , MSE=%.1e' , f1, f2 , zeta1 , zeta2 ,max_isv , max_Fp , mse) 
            % sprintf(' f_1=%.1f , ξ_1=%.2f , i_{sv}= %.1f , MSE=%.1e \n f_2=%.1f , ξ_2=%.2f , F_P=%.1e' , f1,  zeta1 , max_isv ,mse ,f2 ,zeta2 ,max_Fp )
        end
    end
end


%% Plotting response spectre of groung
figure(fig8);
subplot(121)
grid on;
legend();
hold on
semilogx(f_vector, filtered_picos_ddx_ground,'--', 'LineWidth' , 2,'DisplayName', 'Ground');% - Normal

subplot(122)
grid on;
legend();
hold on
semilogx(f_vector, filtered_picos_x_ground,'--', 'LineWidth' , 2,'DisplayName', 'Ground ');%- Normal

%% Save all figures after plotting
% saveas(fig2, 'Input_to_Servo.png');
% saveas(fig7, 'Force_to_Platen.png');
saveas(fig8, 'Response_Spectra.png');


%%

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
    %windowSize = length(f_vector)*0.2; % Adjust window size as needed as percentage of elements in f_vector
    filtered_picos_ddx_m = picos_ddx_m;%movmean(picos_ddx_m, windowSize);
    filtered_picos_x_m  = picos_x_m;%movmean(picos_x_m, windowSize);

end
