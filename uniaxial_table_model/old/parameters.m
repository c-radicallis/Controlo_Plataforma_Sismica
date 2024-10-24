clear all
% Defining parameters
k_p = 1.2993;             % V/cm
tau_sv = 0.0246;          % s
k_sv = 1934.50;        % cm^3/s/V
k_q = 1934.50;        % cm^3/s/V
k_pl = 1.67401e-7;        % m^3/s/kPa  = K_c + C_L
beta_e = 193716.28;       % kPa
V_t = 0.002659;           % m^3
A = 0.012456;             % m^2
m_p = 1.9751;             % ton
m_T = m_p;
c_t = 5.7800;             % kNs/m
K_h = 4*beta_e*A^2/V_t;

%paremeters i couldnt find
% assumed k_c = C_l = k_pl/2
k_c = k_pl/2;
C_l = k_c;

% Converting everything to SI units
k_p  = k_p/100
k_sv =  tau_sv*100^3
k_q = k_q*100^3
k_pl = k_pl/1000


%%%  structure parameters
% 1st mode
m_sp= 2000; % kg
f_sp = 2; % Hz
qsi_sp = 0.05; %

k_sp = m_sp*(2*pi*f_sp)^2; %N/m
c_sp = qsi_sp*2*m_sp*2*pi*f_sp; %N/m/s
%2nd mode
m_sp2= 2000; % kg
f_sp2 = 2; % Hz
qsi_sp2 = 0.05; %

k_sp2 = m_sp*(2*pi*f_sp)^2; %N/m
c_sp2 = qsi_sp*2*m_sp*2*pi*f_sp; %N/m/s


%%

dados = load('elcentro.txt');

t_vector = dados(:,1);
t_step = t_vector(2)
ddx = [t_vector dados(:,2)];
ddy = [t_vector  dados(:,3)];



%%
%{

| **Symbol**   | **Description**                                                     |
|--------------|---------------------------------------------------------------------|
| \( k_p \)    | Proportional gain of the controller                                 |
| \( k_{sv} \) | Servo-valve gain                                                    |
| \( \tau_{sv} \) | Time-delay parameter of the servo-valve transfer function         |
| \( k_q \)    | Valve flow gain                                                     |
| \( k_c \)    | Valve pressure-flow gain                                            |
| \( A \)      | Area of the fluid under compression in the actuator                 |
| \( V_t \)    | Total volume of the fluid under compression in the actuator         |
| \( C_l \)    | Total leakage coefficient of the piston                             |
| \( \beta_e \) | Effective bulk modulus of the system (including oil, entrapped air, etc.) |
| \( k_h \)    | Oil-column frequency of the actuator                                |
| \( m_p \)    | Mass of the platen                                                  |
| \( m_{sp} \) | Mass of the SDOF (Single Degree of Freedom) structure               |
| \( c_t \)    | Combined damping force of the actuator and the platen               |
| \( H_{sp} \) | Transfer function relating the displacements of the platen and the SDOF structure |
| \( k_{sp} \) | Stiffness of the SDOF structure                                     |
| \( c_{sp} \) | Damping of the SDOF structure                                       |
| \( m_t^* \)  | Total mass considering the payload (shaking table)                  |

%}