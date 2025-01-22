clear 

%Controller
k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);

[s,G_T,G_1,G_2,G_T1 ,G_21 ,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT ]=Compute_TFs(G_c);

%%
% sys=G_xT_xref;
% % Convert to zero-pole-gain form
% [z, p, k] = zpkdata(sys, 'v'); % 'v' returns the data in vector format
% zpk_sys = zpk(sys); % Zero-pole-gain form of the system
% 
% % Display the results
% disp('Zeros:');
% disp(z);
% disp('Poles:');
% disp(p);
% disp('Gain:');
% disp(k);
% disp('Zero-Pole-Gain System:');
% disp(zpk_sys);

%%
dados = load('elcentro.txt');
t_vector = dados(:,1);
t_step = t_vector(2)
ddx = [t_vector dados(:,2)];
ddy = [t_vector  dados(:,3)];

%plot(ddx(:,1) ,ddx(:,2))


%% Transfer fuction Bode plots & Output simulation

%  % If axes don't exist, create them
% if ~exist('ax1', 'var') || ~isvalid(ax1)
%     figure(1); % Create a single figure
%     ax1 = subplot(2,2,1); hold(ax1, 'on');
%     ax2 = subplot(2,2,2); hold(ax2, 'on');
%     ax3 = subplot(2,2,3); hold(ax3, 'on');
%     ax4 = subplot(2,2,4); hold(ax4, 'on');
% end
% 
% % First plot
% axes(ax1); % Activate the existing axes
% bode(G_xT_xref);
% title('Bode of G\_xT\_xref'); 
% legend();
% 
% % Second plot
% axes(ax2); % Activate the existing axes
% bode(G_x1_xT);
% title('Bode of G\_x1\_xT'); 
% legend();
% 
% % Third plot
% axes(ax3); % Activate the existing axes
% bode(G_x2_xT);
% title('Bode of G\_x2\_xT'); 
% legend();
% 
% 
% ddx_ref=ddx(:,2);
% ddx_T = lsim(G_xT_xref, ddx_ref , t_vector ,'foh');
% erro = mean((ddx_T-ddx_ref).^2);
% 
% % displacements in milimeters
% x_ref = lsim(1e3/s^2,  ddx(:,2) ,ddx(:,1),'foh');
% x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx(:,2) ,ddx(:,1),'foh');
% 
% 
% % Fourth plot
% axes(ax4); % Activate the existing axes
% title('Acceleration(m/s^2)'); 
% hold on
% plot(ddx(:,1),ddx_ref,"DisplayName","Reference")
% plot(ddx(:,1),ddx_T,"DisplayName","MSE="+string(erro))
% 
% % title('Displacement(mm)'); 
% % plot(ddx(:,1),x_ref,"DisplayName","Reference")
% % plot(ddx(:,1),x_T,"DisplayName","MSE="+string(erro))
% 
% legend()

%% Use PID tuner app to generate a PID controller for the system

% tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
% C = pidtune(G_xT_xref,'PIDF',tuner_opts)
% 
% G_c = C*G_c
% 
% [s,G_T,G_1,G_2,G_T1 ,G_21 ,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT ]=Compute_TFs(G_c)
% 
% % displacements in milimeters
% ddx_T = lsim(G_xT_xref, ddx(:,2) ,ddx(:,1),'foh');
% erro = mean((ddx_T-ddx_ref).^2);
% 
% % x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx(:,2) ,ddx(:,1),'foh');
% % erro = mean((x_T-x_ref).^2);
% 
% axes(ax4); % Activate the existing axes
% hold on
% plot(ddx(:,1),ddx_T,"DisplayName","MSE="+string(erro))
% %plot(ddx(:,1),x_T,"DisplayName","MSE="+string(erro))

%% Finding Response Spectre

f_i=0.1; %freq inicial
f_n=30;  %freq final
f_vector = logspace( log10(f_i) , log10(f_n) , 1000);

m=1;%1kg
zeta=0.05;          %damping ratio (%)

[time_length , n_directions]=size(dados);


picos_ddx_m = zeros( length(f_vector)  , n_directions-1 );
picos_x_m = picos_ddx_m;

for i=1:length(f_vector)    
  
    % from the desired natural frequency, determining stiffness and damping  
    k = m*(2*pi*f_vector(i))^2; %N/m
    c = zeta*2*m*2*pi*f_vector(i); %N/m/s
      
    for j=2:n_directions
        % simulacao do modelo 
        ddx_m = lsim( (c*s+k)/(m*s^2+c*s+k) , dados(:,j) , t_vector ,'zoh');
        picos_ddx_m(i,j-1)=max(abs( ddx_m(:,1) ));  

        x_m = lsim( (c*s+k)/(m*s^2+c*s+k)*1/s^2 , dados(:,j)  , t_vector ,'zoh');
        picos_x_m(i,j-1)=max(abs( x_m(:,1) ));  

    end
end

%valores medios
media_picos_ddx_m=mean(picos_ddx_m,n_directions-1);
media_picos_x_m=mean(picos_x_m,n_directions-1);

%%

% subplot(131)
% plot(dados(:,1),dados(:,2:end))
% title('time series'),xlabel('t(s)'),ylabel('ddx_g (m/s^2)')
% legend()%'ddx' , 'ddy' )% , 'mean')

figure('Position', [100 100 1800 800])
grid on
hold on

subplot(121)
hold on
semilogx(f_vector,picos_ddx_m,'.','MarkerSize',2, 'DisplayName' , 'picos ddx')
semilogx( f_vector, media_picos_ddx_m,'.','MarkerSize',2, 'DisplayName' , 'media picos ddx')
%xlim([0 4])
xlabel('Frequency (Hz)'),ylabel('ddx_m (m/s^2)'),title('Response Spectra')

subplot(122)
hold on
semilogx(f_vector,picos_x_m ,'.','MarkerSize',2, 'DisplayName' , 'picos x')
semilogx(f_vector, media_picos_x_m,'.','MarkerSize',2, 'DisplayName' , 'media picos x')
xlabel('Frequency (Hz)'),ylabel('x_m (m)'),title('Response Spectra')


% Apply a filter to the data

% % Moving average filter
windowSize = 200; % Adjust window size as needed
filtered_picos_ddx_m = movmean(picos_ddx_m, windowSize);
filtered_picos_x_m = movmean(picos_x_m, windowSize);
filtered_media_picos_ddx_m = movmean(media_picos_ddx_m, windowSize);
filtered_media_picos_x_m = movmean(media_picos_x_m, windowSize);

% % Low pass filter
% cut=0.01;
% filtered_picos_ddx_m = lowpass(picos_ddx_m, cut);
% filtered_picos_x_m = lowpass(picos_x_m, cut);

subplot(121)
hold on
semilogx(f_vector, filtered_picos_ddx_m,  'DisplayName' , 'filtrado ddx' , 'LineWidth' , 3)
semilogx(f_vector, filtered_media_picos_ddx_m,  'DisplayName' , 'filtrado media  ddx' , 'LineWidth' , 3)
legend

subplot(122)
hold on
semilogx(f_vector, filtered_picos_x_m,  'DisplayName' , 'filtrado x' , 'LineWidth' , 3)
semilogx(f_vector, filtered_media_picos_x_m,  'DisplayName' , 'filtrado media  x' , 'LineWidth' , 3)
xlim([0 2.5])
legend