function loadTXT(filename)
% loadTXT  Load an LNEC .tgt or .drv text file into MATLAB workspace,
%          creating appropriately named variables in the base workspace.
%
% Usage:
%   loadTXT('LAquilaReducedScale.tgt.txt')
%     → creates variables: time_vector, x_tgt_T, x_tgt_L, x_tgt_V,
%                         ddx_tgt_T, ddx_tgt_L, ddx_tgt_V
%
%   loadTXT('LAquilaReducedScale_34.DRV.txt')
%     → creates variables: x_drv_T_34, x_drv_L_34, x_drv_V_34
%
% The function expects a header block containing at least:
%   Name: <tab‑separated column names>
%   Type: <tab‑separated types ('Displacement' or 'Acceleration')>
%   No. of Samples: <integer>  (only the first value is used)
%   Time step [s]: <float>     (only the first value is used)
%   Data:
%     <tab‑separated numeric rows>
%
% Columns whose data values are all zero will be ignored (not loaded).
%
% For a “.tgt” file, the data columns should contain up to six columns:
%   DispT, DispL, DispV, AccT, AccL, AccV.
% After loading (ignoring zero columns), this creates:
%   time_vector   = (0:(N‑1))' * dt
%   x_tgt_T       = data(:, column where Type=='Displacement' & Name contains 'T')
%   x_tgt_L       = data(:, column where Type=='Displacement' & Name contains 'L')
%   x_tgt_V       = data(:, column where Type=='Displacement' & Name contains 'V')
%   ddx_tgt_T     = data(:, column where Type=='Acceleration' & Name contains 'T')
%   ddx_tgt_L     = data(:, column where Type=='Acceleration' & Name contains 'L')
%   ddx_tgt_V     = data(:, column where Type=='Acceleration' & Name contains 'V')
% Any of these is only created if that column existed and was non‑zero.
%
% For a “.drv” file, the data columns should contain up to three columns:
%   DispT, DispL, DispV.
% After loading (ignoring zero columns), this creates:
%   x_drv_T_i     = data(:, column where Type=='Displacement' & Name contains 'T')
%   x_drv_L_i     = data(:, column where Type=='Displacement' & Name contains 'L')
%   x_drv_V_i     = data(:, column where Type=='Displacement' & Name contains 'V')
% only for those channels that exist (non‑zero).
% “i” is the numeric index parsed from the filename (e.g., “..._34.DRV.txt”).

  %--------------------------%
  % 1. Parse filename parts  %
  %--------------------------%
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s".', ext);
  end

  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    baseName = name;
    suffix   = '';
  else
    baseName = name(1:dotIdx-1);
    suffix   = lower(name(dotIdx+1:end));
  end

  %--------------------%
  % 2. Open the file   %
  %--------------------%
  fid = fopen(filename, 'r');
  if fid < 0
    error('Cannot open "%s" for reading.', filename);
  end

  %----------------------%
  % 3. Read header lines %
  %----------------------%
  numsamps  = [];
  dt        = [];
  colNames  = {};
  typeNames = {};
  while true
    tline = fgetl(fid);
    if ~ischar(tline)
      fclose(fid);
      error('Unexpected end of file while reading header.');
    end
    if startsWith(tline, 'Name:', 'IgnoreCase', true)
      rawNames = strtrim(tline(6:end));
      colNames = strsplit(rawNames, '	');
    elseif startsWith(tline, 'Type:', 'IgnoreCase', true)
      rawTypes = strtrim(tline(6:end));
      typeNames = strsplit(rawTypes, '	');
    elseif startsWith(tline, 'No. of Samples:', 'IgnoreCase', true)
      numsampsArr = sscanf(tline, 'No. of Samples:%f');
      numsamps = numsampsArr(1);
    elseif startsWith(tline, 'Time step [s]:', 'IgnoreCase', true)
      dtArr = sscanf(tline, 'Time step [s]:%f');
      dt = dtArr(1);
    elseif strcmpi(strtrim(tline), 'Data:')
      fgetl(fid);
      break;
    end
  end

  if isempty(numsamps) || isempty(dt) || isempty(colNames) || isempty(typeNames)
    fclose(fid);
    error('Header missing required fields.');
  end

  %------------------------------------%
  % 4. Read numeric data into matrix    %
  %------------------------------------%
  ncol = numel(colNames);
  data = fscanf(fid, repmat('%f',1,ncol), [ncol, Inf])';
  fclose(fid);
  if size(data,1) < numsamps
    warning('Read %d rows but header said %d. Using %d.', size(data,1), numsamps, size(data,1));
    numsamps = size(data,1);
  end

  %--------------------------------------------%
  % 4b. Discard columns that are entirely zero  %
  %--------------------------------------------%
  nzColsMask = any(data ~= 0, 1);
  data       = data(:, nzColsMask);
  colNames   = colNames(nzColsMask);
  typeNames  = typeNames(nzColsMask);

  %--------------------------------------------%
  % 5. Decide behavior based on file suffix    %
  %--------------------------------------------%
  switch suffix
    case 'tgt'
      %---------------------------%
      % Build time vector first  %
      %---------------------------%
      time_vector = (0:(numsamps-1))' * dt;
      assignin('base', 'time_vector', time_vector);

      % Find indices of remaining columns
      dispIdxs = find(strcmpi(typeNames, 'Displacement'));
      accIdxs  = find(strcmpi(typeNames, 'Acceleration'));

      % For each of the six possible channels, attempt to find index
      idxDispT = findChannelIndex(colNames, dispIdxs, 'T');
      idxDispL = findChannelIndex(colNames, dispIdxs, 'L');
      idxDispV = findChannelIndex(colNames, dispIdxs, 'V');
      idxAccT  = findChannelIndex(colNames, accIdxs,  'T');
      idxAccL  = findChannelIndex(colNames, accIdxs,  'L');
      idxAccV  = findChannelIndex(colNames, accIdxs,  'V');

      % Assign only those that were found
      if ~isempty(idxDispT)
        assignin('base', 'x_tgt_T', data(:, idxDispT));
        fprintf('Loaded: x_tgt_T\n');
      end
      if ~isempty(idxDispL)
        assignin('base', 'x_tgt_L', data(:, idxDispL));
        fprintf('Loaded: x_tgt_L\n');
      end
      if ~isempty(idxDispV)
        assignin('base', 'x_tgt_V', data(:, idxDispV));
        fprintf('Loaded: x_tgt_V\n');
      end
      if ~isempty(idxAccT)
        assignin('base', 'ddx_tgt_T', data(:, idxAccT));
        fprintf('Loaded: ddx_tgt_T\n');
      end
      if ~isempty(idxAccL)
        assignin('base', 'ddx_tgt_L', data(:, idxAccL));
        fprintf('Loaded: ddx_tgt_L\n');
      end
      if ~isempty(idxAccV)
        assignin('base', 'ddx_tgt_V', data(:, idxAccV));
        fprintf('Loaded: ddx_tgt_V\n');
      end
      return

    case 'drv'
      % Extract numeric index 'i' from baseName
      tokens = regexp(baseName, '_(\d+)$', 'tokens');
      iStr = tokens{1}{1};

      dispIdxs = find(strcmpi(typeNames, 'Displacement'));
      idxT = findChannelIndex(colNames, dispIdxs, 'T');
      idxL = findChannelIndex(colNames, dispIdxs, 'L');
      idxV = findChannelIndex(colNames, dispIdxs, 'V');

      if ~isempty(idxT)
        var_x_T = sprintf('x_drv_T_%s', iStr);
        assignin('base', var_x_T, data(:, idxT));
        fprintf('Loaded: %s\n', var_x_T);
      end
      if ~isempty(idxL)
        var_x_L = sprintf('x_drv_L_%s', iStr);
        assignin('base', var_x_L, data(:, idxL));
        fprintf('Loaded: %s\n', var_x_L);
      end
      if ~isempty(idxV)
        var_x_V = sprintf('x_drv_V_%s', iStr);
        assignin('base', var_x_V, data(:, idxV));
        fprintf('Loaded: %s\n', var_x_V);
      end
      return

    otherwise
      error('Unsupported file suffix "%s". Expecting tgt or drv (before .txt).', suffix);
  end
end

function idx = findChannelIndex(colNames, idxList, letter)
% findChannelIndex: among idxList, find the first column whose name contains 'letter'
  idx = [];
  for jj = idxList
    if contains(colNames{jj}, letter, 'IgnoreCase', true)
      idx = jj;
      return;
    end
  end
end

