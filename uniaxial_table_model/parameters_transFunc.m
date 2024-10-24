clear all
% Defining parameters

%Servo-valve parameters
tsv=0.0246      %Valve time constant (tsv=0.0246 s)
ksvkq=1934.5    %Valve flow gain (ksvkq=1934.5 cm3/s/V)
kpl=1.67401e-7  %Valve pressure gain & leakadge factor (kpl=1.67401e-7 m3/s/kPa)
Be=193716.28    %Oil Bulk modulus (Be=193716.28 kPa)
Vt=0.002659     %Oil Volume on actuator chamber (Vt=0.002659 m3)
A=0.012456      %Piston area (A=0.012456 m2)
mp=1.9751       %Platen mass (mp=1.9751 t)
ctm=5.78        %Total damping, actuator + platen (ct=5.78 kN s/m)
kh=4*Be*A^2/Vt %(kPa m)

%Additional 1DOF model
msp=4                             %Rigid mass (t)
fn1dof=2                          %Natural frequency (Hz)
zeta1dof=0.02                     %damping ratio (2%)
wn1dof=2*pi*fn1dof           %Natural frequency (rad/s)
cct=2*msp*wn1dof                  %critical damping (t/s)
csp=zeta1dof*cct                  %damping (t/s)
ksp=wn1dof^2*msp                 %stiffness (t/s^2)


%Controller
kp=1.2993       %Pgain (kp=1.2993 V/cm)

s=tf('s')

%1DOF model
Hsp=-msp*s*s/(msp*s*s+csp*s+ksp)
mt=mp+msp*(1+Hsp)
%mt=mp+msp      %massa rigida


Gxr=10^-4*A*kh*kp*ksvkq/(10^-4*A*kh*kp*ksvkq + s*(s*tsv + 1)*(A^2*ctm*s + A^2*kh + A^2*mt*s^2 + kpl*ctm*kh + kpl*kh*mt*s))

%FT do anel aberto Gaa=xp/isv
Gaa=10^-4*A*kh*ksvkq/(s*(s*tsv + 1)*(A^2*ctm*s + A^2*kh + A^2*mt*s^2 + kpl*ctm*kh + kpl*kh*mt*s))

%FT do anel fechado
Gxraf=Gaa*kp/(1+Gaa*kp)


%%

dados = load('elcentro.txt');

t_vector = dados(:,1);
t_step = t_vector(2)
ddx = [t_vector dados(:,2)];
ddy = [t_vector  dados(:,3)];