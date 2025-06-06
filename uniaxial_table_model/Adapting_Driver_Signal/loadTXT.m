function loadTXT(filename)
% loadTXT  Load an LNEC .tgt or .drv text file into MATLAB workspace,
%          creating appropriately named variables in the base workspace.
%
% Usage:
%   loadTXT('LAquilaReducedScale.tgt.txt')
%     → creates variables: time_vector, x_tgt_T, x_tgt_L, x_tgt_V
%
%   loadTXT('LAquilaReducedScale_34.DRV.txt')
%     → creates variables: x_drv_T_34, x_drv_L_34, x_drv_V_34,
%                          ddx_drv_T_34, ddx_drv_L_34, ddx_drv_V_34
%
% The function expects a header block containing at least:
%   Name: <tab‐separated column names>
%   No. of Samples: <integer>  (only the first value is used)
%   Time step [s]: <float>     (only the first value is used)
%   Data:
%     <tab‐separated numeric rows>
%
% For a “.tgt” file, the data columns should contain six columns:
%   DispT, DispL, DispV, AccT, AccL, AccV.
% After loading, this creates:
%   time_vector   = (0:(N–1))' * dt
%   x_tgt_T       = data(:, column labeled "DispT")
%   x_tgt_L       = data(:, column labeled "DispL")
%   x_tgt_V       = data(:, column labeled "DispV")
%
% For a “.drv” file, the data columns should contain six columns:
%   DispT, DispL, DispV, AccT, AccL, AccV.
% (In other words, identical header names, but interpretation differs.)
% After loading, this creates:
%   x_drv_T_i     = data(:, column labeled "DispT")
%   x_drv_L_i     = data(:, column labeled "DispL")
%   x_drv_V_i     = data(:, column labeled "DispV")
%   ddx_drv_T_i   = data(:, column labeled "AccT")
%   ddx_drv_L_i   = data(:, column labeled "AccL")
%   ddx_drv_V_i   = data(:, column labeled "AccV")
% where “i” is the numeric index parsed from the filename (e.g., “…_34.DRV.txt”).
%
% If any expected column is missing, an error is raised.

  %--------------------------%
  % 1. Parse filename parts  %
  %--------------------------%
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s".', ext);
  end

  % Extract the suffix ("tgt" or "DRV") and the baseName (before that)
  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    baseName = name;
    suffix   = '';
  else
    baseName = name(1:dotIdx-1);
    suffix   = lower(name(dotIdx+1:end));  % force lowercase for comparison
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
  numsamps = [];
  dt       = [];
  colNames = {};
  while true
    tline = fgetl(fid);
    if ~ischar(tline)
      fclose(fid);
      error('Unexpected end of file while reading header.');
    end

    if startsWith(tline, 'Name:', 'IgnoreCase', true)
      % Remove "Name:" then trim leading/trailing whitespace
      rawNames = strtrim(tline(6:end));
      % Split on tabs (the file uses tabs between column names)
      colNames = strsplit(rawNames, '\t');
    elseif startsWith(tline, 'No. of Samples:', 'IgnoreCase', true)
      numsampsArr = sscanf(tline, 'No. of Samples:%f'); 
      if isempty(numsampsArr)
        error('Could not parse "No. of Samples" from header.');
      end
      numsamps = numsampsArr(1);  % take the first number
    elseif startsWith(tline, 'Time step [s]:', 'IgnoreCase', true)
      dtArr = sscanf(tline, 'Time step [s]:%f'); 
      if isempty(dtArr)
        error('Could not parse "Time step [s]" from header.');
      end
      dt = dtArr(1);  % take the first value
    elseif strcmpi(strtrim(tline), 'Data:')
      % Next line is repeated "Name" row—skip it
      fgetl(fid);
      break;
    end
  end

  if isempty(numsamps) || isempty(dt) || isempty(colNames)
    fclose(fid);
    error('Header did not contain required fields (Name, No. of Samples, Time step).');
  end

  %-----------------------------------%
  % 4. Read numeric data into matrix  %
  %-----------------------------------%
  ncol = numel(colNames);
  fmt  = repmat('%f', 1, ncol);
  data = fscanf(fid, fmt, [ncol, Inf])';
  fclose(fid);

  % If fewer rows were read than expected, adjust numsamps
  if size(data,1) < numsamps
    warning('Read %d rows but header said %d. Using %d.', size(data,1), numsamps, size(data,1));
    numsamps = size(data,1);
  end

  %--------------------------------------------%
  % 5. Decide behavior based on file suffix    %
  %--------------------------------------------%

  switch suffix
    case 'tgt'
      %--------------------------%
      % 5a. .tgt  -> Displacements only  %
      %--------------------------%
      % Look up exact column names "DispT", "DispL", "DispV"
      idxT = find(strcmpi(colNames, 'DispT'), 1);
      idxL = find(strcmpi(colNames, 'DispL'), 1);
      idxV = find(strcmpi(colNames, 'DispV'), 1);

      if isempty(idxT) || isempty(idxL) || isempty(idxV)
        error('Could not find DispT, DispL, or DispV columns in .tgt file header.');
      end

      % Build time vector
      time_vector = (0:(numsamps-1))' * dt;

      % Assign to base workspace
      assignin('base', 'time_vector', time_vector);
      assignin('base', 'x_tgt_T',       data(:, idxT));
      assignin('base', 'x_tgt_L',       data(:, idxL));
      assignin('base', 'x_tgt_V',       data(:, idxV));

      fprintf('Loaded variables into base workspace:\n');
      fprintf('  time_vector\n');
      fprintf('  x_tgt_T\n');
      fprintf('  x_tgt_L\n');
      fprintf('  x_tgt_V\n');
      return

    case 'drv'
      %--------------------------------%
      % 5b. .drv -> Displacements & Accelerations  %
      %--------------------------------%
      % Extract numeric index 'i' from baseName (expects “…_<i>”)
      tokens = regexp(baseName, '_(\d+)$', 'tokens');
      if isempty(tokens)
        error('Could not parse numeric index from baseName "%s".', baseName);
      end
      iStr = tokens{1}{1};  % e.g. '34'

      % Find exact header labels
      idxDispT = find(strcmpi(colNames, 'DispT'), 1);
      idxDispL = find(strcmpi(colNames, 'DispL'), 1);
      idxDispV = find(strcmpi(colNames, 'DispV'), 1);
      idxAccT  = find(strcmpi(colNames, 'AccT'),  1);
      idxAccL  = find(strcmpi(colNames, 'AccL'),  1);
      idxAccV  = find(strcmpi(colNames, 'AccV'),  1);

      if any(cellfun(@isempty, {idxDispT, idxDispL, idxDispV, idxAccT, idxAccL, idxAccV}))
        error(['Could not find one or more of DispT/DispL/DispV/AccT/AccL/AccV ', ...
               'columns in .drv file header.']);
      end

      % Build variable names with index suffix
      var_x_T   = sprintf('x_drv_T_%s',   iStr);
      var_x_L   = sprintf('x_drv_L_%s',   iStr);
      var_x_V   = sprintf('x_drv_V_%s',   iStr);
      var_ddx_T = sprintf('ddx_drv_T_%s', iStr);
      var_ddx_L = sprintf('ddx_drv_L_%s', iStr);
      var_ddx_V = sprintf('ddx_drv_V_%s', iStr);

      % Assign into base workspace
      assignin('base', var_x_T,   data(:, idxDispT));
      assignin('base', var_x_L,   data(:, idxDispL));
      assignin('base', var_x_V,   data(:, idxDispV));
      assignin('base', var_ddx_T, data(:, idxAccT));
      assignin('base', var_ddx_L, data(:, idxAccL));
      assignin('base', var_ddx_V, data(:, idxAccV));

      fprintf('Loaded variables into base workspace:\n');
      fprintf('  %s\n', var_x_T);
      fprintf('  %s\n', var_x_L);
      fprintf('  %s\n', var_x_V);
      fprintf('  %s\n', var_ddx_T);
      fprintf('  %s\n', var_ddx_L);
      fprintf('  %s\n', var_ddx_V);
      return

    otherwise
      error('Unsupported file suffix "%s". Expecting .tgt or .drv (before .txt).', suffix);
  end
end
