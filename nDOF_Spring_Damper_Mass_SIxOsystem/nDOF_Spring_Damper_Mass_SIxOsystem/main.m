clc; clear all; 
%%
n= 2;
M= ones(1,n)*50;
C= ones(1,n)*10;
K= ones(1,n)*20000;

%% input signal
in_opt.f= 5; % frequency in Hz
in_opt.a0= 100; % force amplitude in N
%% time 
t= 0:0.01:1;

%% calculating state space matrices A and B
[A,B]= Model_nDOF_Spring_Mass_Damper_SIxOsystem_Force(n,M,C,K);

%% calculating state space matrix C
% considering the position of the last mass as output
C_=[zeros(1,n-1) 1 zeros(1,n)];
% considering the position of the first mass as output
% C_=[1 zeros(1,n-1) zeros(1,n)];
%% calculating state space matrix D
D= zeros;

input_signal= @(t) in_opt.a0*sin(2*pi*in_opt.f*t);
fun= @(t,x) A*x+B* input_signal(t);

x0= zeros(size(A,1),1); % initial values
[t_out,x_out]= ode45(fun,t,x0);

y= C_*x_out';

plot(t,y)
[~,body_considered]= max(C_);
title(['Position of body ' num2str(body_considered) ' vs Time' ])
xlabel('Time (sec)')
ylabel('Position (m)')

