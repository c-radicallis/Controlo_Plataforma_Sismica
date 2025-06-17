% Filename: plot_FRF_bode.m
% Description:
%   Reads lines 2–14 from 'identified_FRF.txt', extracts column 1 (frequency [Hz])
%   and column 2 (magnitude, linear), then plots:
%     - Bode plot (FRD with zero-phase)
%     - Magnitude-only semilogx plot (20*log10(mag) vs freq)
%
% Usage:
%   Place this file in your MATLAB path or current folder. Make sure
%   'identified_FRF.txt' is in the current folder (or adjust filename).
%   Then in MATLAB: >> plot_FRF_bode
%
% Note:
%   - Adjust `filename`, number of data lines (currently 13 rows after header),
%     or delimiter if your file uses different separators.
%   - If your data file has more or fewer rows or a different structure, update
%     the section marked “%% Adjust if needed” accordingly.


%% 0. Clear and define filename
clear;
clc;
close all;

filename = 'identified_FRF.txt';  % <-- change if your file has a different name or path

%% 1. Open file and skip header
fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open file "%s". Make sure the file exists in the current folder or provide a correct path.', filename);
end

% Read and discard the first line (header)
headerLine = fgetl(fid);
if ~ischar(headerLine)
    fclose(fid);
    error('File "%s" appears empty or unreadable.', filename);
end

%% 2. Read the next 13 lines: two numeric columns, skip rest of each line
%    Adjust the count 13 if you have a different number of data rows.
numDataRows = 13;  % lines 2 through 14 inclusive

% The format '%f%f%*[^\n]' reads:
%   %f : first numeric field (frequency)
%   %f : second numeric field (magnitude)
%   %*[^\n] : skip the rest of the line until newline
% By default textscan uses whitespace (spaces/tabs) as delimiters for %f.
% If your file strictly uses tabs and you want to ensure that, you can add:
%   'Delimiter','\t'
% inside textscan. Usually '%f' handles whitespace-separated numbers fine.
data = textscan(fid, '%f%f%*[^\n]', numDataRows);
fclose(fid);

% Extract vectors
freq = data{1};       % frequency in Hz, expected size [numDataRows×1]
magLinear = data{2};  % linear magnitude, same size

%% 3. Validate read results
if numel(freq) < numDataRows || numel(magLinear) < numDataRows
    warning('Only %d rows read, expected %d. Check file formatting or adjust numDataRows.', numel(freq), numDataRows);
end

% Check for NaNs in parsed data
nanIdx = isnan(freq) | isnan(magLinear);
if any(nanIdx)
    warning('Found NaN in parsed data at rows: %s. Inspect your file for missing or non-numeric entries.', mat2str(find(nanIdx)'));
    % Display problematic rows
    disp(table(find(nanIdx), freq(nanIdx), magLinear(nanIdx), ...
         'VariableNames', {'RowIndex','freq','magLinear'}));
    % You may choose to remove or interpolate NaNs here. For now, we stop if NaNs exist:
    error('Cannot proceed with NaN entries. Please correct or remove invalid data lines.');
end

%% 4. Build FRD model assuming zero phase, then plot Bode
% Convert frequency from Hz to angular frequency (rad/s)
w = 2*pi * freq;

% Build complex frequency response H(jw) with zero phase
H = magLinear .* exp(1j * zeros(size(freq)));

% Create FRD model
sys = frd(H, w);

% Plot Bode
figure('Name','Bode Plot (zero-phase assumption)','NumberTitle','off');
bode(sys);
grid on; 
title('Bode Plot (Magnitude from file, Phase assumed zero)');

% %% 5. Plot magnitude-only semilogx (20*log10) vs frequency
% figure('Name','Magnitude vs Frequency','NumberTitle','off');
% semilogx(freq, 20*log10(magLinear), '-o', 'LineWidth', 1.5);
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Magnitude vs Frequency (Lines 2–14 of data file)');

%% 6. End of script
% You can save or export figures as needed.
% If you obtain phase data later, you can replace the zero-phase assumption by:
%   H = magLinear .* exp(1j * deg2rad(phaseDeg));
% where phaseDeg is a numeric vector of same size.

%% ================================================
% Alternative approach using readtable (commented out)
% If you prefer to use readtable and your file is tab-delimited with headers:
% Uncomment and adjust variable names as needed.

%{
% 1) Detect import options
opts = detectImportOptions(filename, 'FileType','text','Delimiter','\t');
opts.VariableNamesLine = 1;   % header row
opts.DataLines = [2, 14];     % lines 2 through 14 inclusive

% 2) Inspect sanitized names
disp('Sanitized variable names:');
disp(opts.VariableNames);
disp('Original headers (descriptions):');
disp(opts.VariableDescriptions);

% Suppose sanitized names are:
%   col1 = opts.VariableNames{1};  % e.g. 'Frequency__Hz___Plot_0'
%   col2 = opts.VariableNames{2};  % e.g. 'Magnitude___Plot_0'
col1 = opts.VariableNames{1};
col2 = opts.VariableNames{2};

% 3) Force these two to double and select only them
opts = setvartype(opts, {col1, col2}, 'double');
opts.SelectedVariableNames = {col1, col2};

% 4) Read table
T = readtable(filename, opts);

% 5) Extract numeric arrays
freq_tbl = T.(col1);
magLinear_tbl = T.(col2);

% 6) Check for NaNs
if any(isnan(freq_tbl)) || any(isnan(magLinear_tbl))
    error('NaN entries detected in table import. Fix data or import options.');
end

% 7) Plot as above
w_tbl = 2*pi * freq_tbl;
H_tbl = magLinear_tbl .* exp(1j * zeros(size(freq_tbl)));
sys_tbl = frd(H_tbl, w_tbl);
figure; bode(sys_tbl); grid on; title('Bode via readtable');
figure; semilogx(freq_tbl, 20*log10(magLinear_tbl), '-o'); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude vs Frequency via readtable');
%}

% End of alternative approach

%% Loading Model with Standard Tune

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

% s,G_T,G_1,G_2,G_T1 ,G_21 ,G_svq,G_csv,G_x2_x1,G_x1_xT,G_xT_Fp,G_Fp_xref,G_xT_xref,G_x1_xref,G_x2_xT , G_Fp_isv  ,c1,c2,k1,k2,AA , BB , CC , DD 
[~,~,~,~,~ ,~ ,~,~,~,~,~,~,G_xT_xref,~,~ , ~ ,~,~,~,~ ,~ , ~ , ~ , ~  ]=Compute_TFs(G_c, mT , cT , m1 , m2 , f1, zeta1 , f2 , zeta2);

hold on;
bode(G_xT_xref)

%%

%% 0. Clear and define filename
clear;
clc;
close all;

filename = 'identified_FRF.txt';  % <-- change if your file has a different name or path

headerLines = 1;        % number of header lines to skip (your first line is the header)

% --- Step 1: Read the data ---
% Option A: using readmatrix (R2020b or later)
% This will read numeric data; blank/missing entries become NaN.
dataAll = readmatrix(filename, 'FileType', 'text', 'Delimiter', '\t', 'NumHeaderLines', headerLines);

% dataAll may have more columns (e.g., 4 columns in header but only first two numeric). 
% Extract the first two columns which correspond to Frequency [Hz] and Magnitude.
% Remove rows that are entirely NaN or where first column is NaN:
validRows = ~isnan(dataAll(:,1));  
data = dataAll(validRows, 1:2);

% Option B: using readtable + detectImportOptions (more control over variable names)
% opts = detectImportOptions(filename, 'Delimiter', '\t');
% opts.DataLines = [2 Inf];   % start reading from line 2 onward
% % If the file has extra blank columns you don’t need, select only first two variables:
% if numel(opts.VariableNames) >= 2
%     opts.SelectedVariableNames = opts.VariableNames(1:2);
% end
% T = readtable(filename, opts);
% data = T{:,1:2};  % numeric array: col1 = frequency, col2 = magnitude

% Now data is an N×2 array. If you only want lines 2–14 specifically and file may contain more:
% Suppose dataAll includes exactly lines 2..end; then:
% nRows = size(data,1);
% If nRows >= 13, keep only first 13 rows:
if size(data,1) >= 13
    data = data(1:13, :);
end

% --- Step 2: Separate frequency and magnitude ---
freq = data(:,1);   % in Hz
mag_lin = data(:,2);  % linear magnitude

% --- (Optional) Inspect raw values ---
% disp(table(freq, mag_lin));

% --- Step 3: Convert magnitude to dB for Bode plot ---
% Avoid log of zero or negative: ensure mag_lin > 0
if any(mag_lin <= 0)
    warning('Some magnitude values are <= 0; these cannot be converted to dB. They will be set to -Inf or NaN.');
end
mag_dB = 20*log10(mag_lin);

% --- Step 4: Plotting ---
figure;
semilogx(freq, mag_dB, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); 
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Bode‐style Magnitude Plot (data from lines 2–14)');
% Optionally adjust axis limits:
xlim([min(freq)*0.8, max(freq)*1.2]);  
% ylim can be set if desired, e.g.: ylim([min(mag_dB)-5, max(mag_dB)+5]);

% If you also want to plot the linear magnitude (on linear y-scale) versus log-x:
figure;
semilogx(freq, mag_lin, 's-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (linear)');
title('Magnitude vs Frequency (linear scale)');

% --- Notes ---
% 1. If your file uses spaces instead of tabs, change 'Delimiter','\t' accordingly, e.g. readmatrix will auto-detect spaces.
% 2. If you only want lines 2 to 14 exactly, ensure that readmatrix/readtable correctly skips header and then you index rows 1:13 of the imported numeric block.
% 3. If you prefer plain “plot” on a linear x-axis (not typical for Bode), replace semilogx with plot.
% 4. If you have phase data in another column, you could similarly read it and plot with semilogx, but here only magnitude is shown.
% 5. Ensure your file encoding and delimiters match: you can open it in a text editor to confirm tabs/spaces. You can also use “importdata” or “textscan” if you need finer control:
%
%    % Example with textscan:
%    fid = fopen(filename, 'r');
%    fgetl(fid);  % skip header line
%    C = textscan(fid, '%f%f%*s%*s', 'Delimiter', '\t'); 
%      % '%f%f%*s%*s' reads two floats, ignores next two fields per row if present
%    fclose(fid);
%    freq = C{1};
%    mag_lin = C{2};
%    if length(freq) >= 13
%        freq = freq(1:13);
%        mag_lin = mag_lin(1:13);
%    end
%    mag_dB = 20*log10(mag_lin);
%    figure; semilogx(freq, mag_dB, 'o-'); grid on;
%
% 6. Place the MATLAB script (or function) in the same folder as the data file or supply full path.
% 7. If you run into NaNs, display dataAll in the workspace and check which columns contain numeric entries.
