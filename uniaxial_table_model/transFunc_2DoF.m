%clear all

%Controller
k_p=1.2993/1e-2 %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1)  

[s,G_T,G_1,G_2,G_T1 ,G_21 ,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT ]=Compute_TFs(G_c)

%%
dados = load('elcentro.txt');

t_vector = dados(:,1);
t_step = t_vector(2)
ddx = [t_vector dados(:,2)];
ddy = [t_vector  dados(:,3)];

%plot(ddx(:,1) ,ddx(:,2))


%%
% If axes don't exist, create them
if ~exist('ax1', 'var') || ~isvalid(ax1)
    figure(1); % Create a single figure
    ax1 = subplot(2,2,1); hold(ax1, 'on');
    ax2 = subplot(2,2,2); hold(ax2, 'on');
    ax3 = subplot(2,2,3); hold(ax3, 'on');
    ax4 = subplot(2,2,4); hold(ax4, 'on');
end

% First plot
axes(ax1); % Activate the existing axes
bode(G_xT_xref);
title('Bode of G\_xT\_xref'); 
legend();

% Second plot
axes(ax2); % Activate the existing axes
bode(G_x1_xT);
title('Bode of G\_x1\_xT'); 
legend();

% Third plot
axes(ax3); % Activate the existing axes
bode(G_x2_xT);
title('Bode of G\_x2\_xT'); 
legend();


ddx_ref=ddx(:,2);
ddx_T = lsim(G_xT_xref, ddx(:,2) ,ddx(:,1),'foh');
erro = mean((ddx_T-ddx_ref).^2);

% displacements in milimeters
x_ref = lsim(1e3/s^2,  ddx(:,2) ,ddx(:,1),'foh');
x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx(:,2) ,ddx(:,1),'foh');


% Fourth plot
axes(ax4); % Activate the existing axes
title('Acceleration(m/s^2)'); 
hold on
plot(ddx(:,1),ddx_ref,"DisplayName","Reference")
plot(ddx(:,1),ddx_T,"DisplayName","MSE="+string(erro))

% title('Displacement(mm)'); 
% plot(ddx(:,1),x_ref,"DisplayName","Reference")
% plot(ddx(:,1),x_T,"DisplayName","MSE="+string(erro))

legend()




%%
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
C = pidtune(G_xT_xref,'PIDF',tuner_opts)

G_c = C*G_c

[s,G_T,G_1,G_2,G_T1 ,G_21 ,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT ]=Compute_TFs(G_c)

% displacements in milimeters
ddx_T = lsim(G_xT_xref, ddx(:,2) ,ddx(:,1),'foh');
erro = mean((ddx_T-ddx_ref).^2);

% x_T = lsim(G_xT_xref*1e3/s^2 ,  ddx(:,2) ,ddx(:,1),'foh');
% erro = mean((x_T-x_ref).^2);

axes(ax4); % Activate the existing axes
hold on
plot(ddx(:,1),ddx_T,"DisplayName","MSE="+string(erro))
%plot(ddx(:,1),x_T,"DisplayName","MSE="+string(erro))
