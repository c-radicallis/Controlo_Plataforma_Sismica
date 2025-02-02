clear
clc

%coupled 2DOF system
%  structure parameters
% 1st mode
m1= 1e3; % kg
f1 = 2; % Hz
zeta1 = 0.02 ; 
% k1 = m1*(2*pi*f1)^2; %N/m
c1 = zeta1*2*m1*2*pi*f1; %N/m/s

%2nd mode
m2= 1e3; % kg
f2 = 10; % Hz
zeta2 = 0.05; %
c2 = zeta2*2*m2*2*pi*f2; %N/m/s

syms k1 k2
assume(k1 ,"positive")
assume(k2 ,"positive")

Y = vpasolve([
    (m1*k2 + m2*(k1+k2))/(2*m1*m2) - 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f1)^2  ,
    (m1*k2 + m2*(k1+k2))/(2*m1*m2)+ 0.5*sqrt( ((m1*k2+m2*(k1+k2))/(m1*m2))^2 - 4*( k1*k2 )/(m1*m2) ) == (2*pi*f2)^2,
], [k1,k2])

k1 = double(Y.k1)
k2 = double(Y.k2)

% K=[sol_k1+sol_k2 -sol_k2;-sol_k2 sol_k2];
% M=diag([m1,m2]);
% [Vectors , Values] = eig(K/M)
% 
% f_1=double(sqrt(Values(1,1))/(2*pi))
% f_2=double(sqrt(Values(2,2))/(2*pi))

