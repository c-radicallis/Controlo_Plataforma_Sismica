clc
clear
close all

s = tf('s');

%% Systems
% Nominal system
G_0 = 0.1/(s^2 + 1.5*s + 1);
% Perturbed system
G_err = 0.9/(s^2 + 1.5*s + 1);

% Controllers
C_OL = 10;
C_CL = 55 + 30/s + 25*s;

% Open-loop systems
OL_0    = C_OL*G_0;
OL_Gerr = C_OL*G_err;

% Closed-loop transfer functions (reference -> output)
T_0    = feedback(C_CL*OL_0,1);
T_Gerr = feedback(C_CL*OL_Gerr,1);


%% Time vector
t = 0:0.01:20; t=t';

%% Reference signal
% r(t=0)=0
% r(0<t<5)=5
% r(5<t<10)=10
% r(t>10)=0
r = zeros(size(t));

r(t > 0  & t < 5)  = 5;
r(t >= 5 & t < 10) = 10;
r(t >= 10)         = 0;

%% System responses

% Output responses
y0_CL   = lsim(T_0, r, t);
yerr_CL = lsim(T_Gerr, r, t);
y0_dist = y0 + 0.5*sin(2*pi*0.5*t);% Disturbed output system

% Control actions
u0_OL   = C_OL*y0;

u0_CL=lsim(C_CL,r-y0 , t);
uerr_CL=lsim(C_CL,r-yerr , t);

% % Closed-loop transfer functions (reference -> control action)
% U_0    = feedback(C_CL, G_0);
% U_Gerr = feedback(C_CL, G_err);
% % Control actions
% u0   = lsim(U_0, r, t);
% uerr = lsim(U_Gerr, r, t);

%% Top plot: reference and outputs
figure
subplot(2,1,1)

plot(t,r,'k--','LineWidth',1.5)
hold on

plot(t,y0,'b','LineWidth',2)
plot(t,yerr,'r','LineWidth',2)
plot(t,y0_dist,'g','LineWidth',1.5)

grid on

xlabel('Time (s)')
ylabel('Output')

title('Reference and System Outputs')

legend('Reference', ...
       'Nominal System', ...
       'Perturbed System', ...
       'Nominal + Sinusoidal Disturbance')

%% Bottom plot: control actions

subplot(2,1,2)

plot(t,u0,'b','LineWidth',2)
hold on

plot(t,uerr,'r','LineWidth',2)

grid on

xlabel('Time (s)')
ylabel('Control Action')

title('Control Inputs')

legend('Nominal System', ...
       'Perturbed System')