clear;clc;close all;

addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'\Adapting_Driver_Signal\PRJ_project\
%% Load target

% Converting the .drv to .txt, and loading the txt
script_folder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\personal_python_packages';  % modify to your actual path
if count(py.sys.path, script_folder) == 0
    insert(py.sys.path, int32(0), script_folder);
end

base_name = 'LAquilaReducedScale_0';   % or get from user input % 2. Define only the name (no folder); you can prompt the user or set it manually:
ext = '.drv';                         % driver extension

input_folder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project'; % 1. Define folder once
in_file = fullfile(input_folder, [base_name, ext]); % 3. Create full path
if ~isfile(in_file) % 4. (Optional) Check existence before proceeding
    error('Input file does not exist: %s', in_file);
end

out_dir = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project';

try
    py_output = py.LTF_to_TXT.ltf_to_txt(in_file, out_dir);     % py_output is a Python string; convert to MATLAB char:
    output_path = char(py_output);
    fprintf('Python function returned output path: %s\n', output_path);
catch ME
    disp('Error calling Python function:');
    disp(ME.message);
end

loadTXT('LAquilaReducedScale_0.drv.txt')


%% Finding Target Response Spectre
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_tgt , picos_x_tgt] = ResponseSpectrum( time_vector , x_tgt , ddx_tgt, f_vector , 1);

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal


%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

mT=1.9751*1e3; %Platen mass (mp=1.9751 t)
cT=5.78*1e3;   %Total damping, actuator + platen (ct=5.78 kN s/m1)
mass=2e3;
m1 = mass; % kg % 1st mode
f1 = 1.5; % Hz   % 1.5 < f1 < 4
zeta1 = 0.02 ; % 2 < zeta1 < 10
m2 = mass; % kg %2nd mode
f2 =6; % Hz % 6 < f2 < 10
zeta2 = 0.05; % 5 < zeta2 < 25r

c1 = zeta1*2*m1*2*pi*f1; c2 = zeta2*2*m2*2*pi*f2; %N/m/s%coupled 2DOF system

k_p=1.2993/1e-2; %SI units %Pgain (kp=1.2993 V/cm) 
G_c = tf(k_p,1);% Controller

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2, ss_model 
[~,~,~,~,~ ,~ ,~,~,~,~,~,~,G_xT_xref,~,~ , ~ ,~,~,~,~ , ~ , ~ , ~ , ~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

%% Finding Response Spectre of Ground
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_tgt , picos_x_tgt] = ResponseSpectrum( t_vector , x_tgt , ddx_tgt, f_vector , 1);

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt(:, 1),'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
