clear;clc;close all;


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


%% 

% Suppose your ltf_utils.py is in '/path/to/py-scripts'
script_folder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\personal_python_packages';  % modify to your actual path

if count(py.sys.path, script_folder) == 0
    insert(py.sys.path, int32(0), script_folder);
end

% Paths as strings; MATLAB will convert to Python str automatically.
in_file ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\LNEC_Adapta_Driver\LNEC_ERIES_RE-SAFE\CTL\SystemId\pink_noise_40Hz_T3mm.drv';
out_dir = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project';

try
    py_output = py.LTF_to_TXT.ltf_to_txt(in_file, out_dir);
    % py_output is a Python string; convert to MATLAB char:
    output_path = char(py_output);
    fprintf('Python function returned output path: %s\n', output_path);
catch ME
    disp('Error calling Python function:');
    disp(ME.message);
end

%%

oldName = fullfile(out_dir , 'pink_noise_40Hz_T3mm.drv.txt');
newName = fullfile(out_dir , 'pink_noise_40Hz_T3mm_0.drv.txt');
movefile(oldName, newName)


%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project\'
filename = 'pink_noise_40Hz_T3mm_0.drv.txt';
loadTXT(filename)

%%   % --- Simulation --
t_vector = time_drv_0;
t_step = t_vector(2);
x_drv    = x_drv_T_0;
ddx_drv = secondDerivativeTime(x_drv,t_step);

x_acq = lsim(G_xT_xref ,  x_drv ,t_vector,'foh');
ddx_acq = secondDerivativeTime(x_acq,t_step);

 hold on; grid on; legend()
 plot(t_vector , x_drv)
 plot(t_vector,x_acq)

%%  --- Write .txt file from simulated "acquired" data 
filename_acq = strrep(filename, '.drv.txt', '_acq.txt');
save_folder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_project'
 writeTXT(t_vector , x_acq , ddx_acq , save_folder , filename_acq)

%% Run .txt to .LTF conversion

full_filename_acq = fullfile(save_folder, filename_acq);
try
    py_output = py.TXT_to_LTF.txt_to_ltf(full_filename_acq, out_dir);
    % py_output is a Python string; convert to MATLAB char:
    output_path = char(py_output);
    fprintf('Python function returned output path: %s\n', output_path);
catch ME
    disp('Error calling Python function:');
    disp(ME.message);
end
