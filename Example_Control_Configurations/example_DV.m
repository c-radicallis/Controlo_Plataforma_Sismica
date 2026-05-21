% Example 10.7.5 - from Van den Hof lecture notes
% two stage closed loop identification method modified for
% reference r2 only, and proportional controller

clear;clc; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-leftwin_size
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
close all;

s = tf('s'); Ts=0.01;

%% Time vector
t = [0:Ts:4]';
v = 1*sin(2*pi*1*t);

r2 = zeros(size(t));
r2(t > 0  & t < 2)  = 5;
r2(t >= 2 & t < 4) = 0;
r2(t >= 20)         = 0;

%% Systems
% Nominal system
G_0 = (s+1)/(s*(s+10));%0.1/(s^2 + 1.5*s + 1);
% Model error system
G_err = 1.1*G_0

% Controllers
% C_OL = tf(10);
C_CL = pid(10,50,50,1)
C_OL = 10*s/(s+1);

% Open-loop systems
OL_0    = C_OL*G_0;
OL_Gerr = C_OL*G_err;

%% Open loop sim
y_0_OL=lsim(OL_0,r2,t);
y_err_OL = lsim(OL_Gerr, r2, t);
y_0_dist_OL = y_0_OL+v;

u_0_OL=lsim(C_OL,r2,t);

%% Closed loop sim
r =  lsim( C_CL , r2, t );%reference increases in magnitude when controller gain decreases
% default
S_0 = 1/(1+G_0*C_CL);
y_0_CL=lsim(G_0*S_0 , r , t);
u_0_CL = r - lsim(C_CL , y_0_CL ,t);

%model error
S_err = 1/(1+G_err*C_CL);
y_err_CL=lsim(G_err*S_err , r , t);
u_err_CL = r - lsim(C_CL , y_err_CL ,t);

% disturbance
y_0_dist_CL=lsim(G_0*S_0 , r , t)+lsim(S_0 , v , t);
u_0_dist_CL = r - lsim(C_CL , y_0_dist_CL ,t);


% combined control
% C_OL = 10*s/(s+1);
% C_CL = pid(20,20,20,1)
r1 = lsim(C_OL , r2 ,t);
r = r1 + lsim( C_CL , r2 , t );

% default
y_0_Comb=lsim(G_0*S_0 , r , t);
u_0_Comb = r - lsim(C_CL , y_0_Comb ,t);

%model error
y_err_Comb=lsim(G_err*S_err , r , t);
u_err_Comb = r - lsim(C_CL , y_err_Comb ,t);

% disturbance
y_0_dist_Comb=lsim(G_0*S_0 , r , t)+lsim(S_0 , v , t);
u_0_dist_Comb = r - lsim(C_CL , y_0_dist_Comb ,t);



%% Plotting
% MATLAB default color palette
colors = get(groot,'defaultAxesColorOrder');
c_nom  = colors(1,:);   % blue
c_err  = colors(2,:);   % orange
c_dist = colors(3,:);   % yellow
% blend = 0; 
% c_nom  =  c_nom + blend*(1-c_nom);
% c_err =  c_err + blend*(1-c_err);
% c_dist = c_dist + blend*(1-c_dist);
%% Top plot - Outpu
fig1=figure(1);
subplot(2,1,1)
plot(t,r2,'k-','LineWidth',1.8)
hold on

% Open-loop responses
plot(t,y_0_OL     ,'Color',c_nom ,'LineWidth',3)
plot(t,y_err_OL  ,'-.' ,'Color',c_err ,'LineWidth',3)
plot(t,y_0_dist_OL,'Color',c_dist,'LineWidth',3)
ylim([-1,6])
grid on

xlabel('Time (s)')
ylabel('Output')

title('Reference and System Outputs')

legend( ...
    'Reference', ...
    'OL: Nominal', ...
    'OL: Model error', ...
    'OL: Disturbance', ...
    'Location','best')

%% Bottom plot - Control actions

subplot(2,1,2)
% Open-loop control
plot(t,u_0_OL,'-','Color',c_nom,'LineWidth',3);
hold on
ylim([-200,200])

grid on

xlabel('Time (s)')
ylabel('Control Action')

title('Control Inputs')

fontsize(scale=2); set(fig1,'WindowState', 'maximized');

%% Top plot - Outpu
fig2=figure(2);
subplot(2,1,1)
plot(t,r2,'k-','LineWidth',1.8)
hold on

% Closed-loop responses
plot(t,y_0_CL        ,'-','Color',c_nom   ,'LineWidth',3)
plot(t,y_err_CL     ,'-.','Color',  c_err ,'LineWidth',3)
plot(t,y_0_dist_CL,'-','Color',c_dist ,'LineWidth',3)
ylim([-1,6])
grid on

xlabel('Time (s)')
ylabel('Output')

title('Reference and System Outputs')

legend( ...
    'Reference', ...
    'CL: Nominal', ...
    'CL: Model error', ...
    'CL: Disturbance', ...
    'Location','best')

%% Bottom plot - Control actions

subplot(2,1,2);hold on
% Closed-loop control
plot(t,u_0_CL     ,'-.',     'Color',c_nom ,'LineWidth',3)
plot(t,u_err_CL   ,'-.',   'Color',c_err ,'LineWidth',3)
plot(t,u_0_dist_CL,'-','Color',c_dist  ,'LineWidth',3)
ylim([-200,200])

grid on

xlabel('Time (s)')
ylabel('Control Action')
title('Control Inputs')

legend( ...
    'CL: Nominal', ...
    'CL: Model error', ...
    'CL: Disturbance', ...
    'Location','best')

fontsize(scale=2); set(fig2,'WindowState', 'maximized');

%%op Combined control plot - Outpu
fig3=figure(3);
subplot(2,1,1)
plot(t,r2,'k-','LineWidth',1.8)
hold on

% Closed-loop responses
plot(t,y_0_Comb        ,'-','Color',c_nom   ,'LineWidth',3)
plot(t,y_err_Comb     ,'-.','Color',  c_err ,'LineWidth',3)
plot(t,y_0_dist_Comb,'-','Color',c_dist ,'LineWidth',3)
ylim([-1,6])
grid on

xlabel('Time (s)')
ylabel('Output')

title('Reference and System Outputs')

legend( ...
    'Reference', ...
    'Combined: Nominal', ...
    'Combined: Model error', ...
    'Combined: Disturbance', ...
    'Location','best')

%Bottom plot - Control actions

subplot(2,1,2);hold on
plot(t,u_0_Comb     ,'-',     'Color',c_nom ,'LineWidth',3)
plot(t,u_err_Comb  ,'-.',   'Color',c_err ,'LineWidth',3)
plot(t,u_0_dist_Comb,'-','Color',c_dist  ,'LineWidth',3)
ylim([-200,200])
grid on

xlabel('Time (s)')
ylabel('Control Action')
title('Control Inputs')

legend( ...
    'Combined: Nominal', ...
    'Combined: Model error', ...
    'Combined: Disturbance', ...
    'Location','best')

fontsize(scale=2); set(fig3,'WindowState', 'maximized');