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
%   No. of Samples: <integer>
%   Time step [s]: <floating point>
%   Data:
%     <tab‐separated numeric rows>
%
% For a “.tgt” file, the data columns should contain exactly the
% three displacement channels T, L, V (possibly with zero columns
% that will be discarded).  After loading, this creates:
%   time_vector   = (0:(N–1))' * dt
%   x_tgt_T       = column corresponding to T displacement
%   x_tgt_L       = column corresponding to L displacement
%   x_tgt_V       = column corresponding to V displacement
%
% For a “.drv” file, the data columns should contain six non‐zero
% channels: three displacements (T, L, V) and three accelerations
% (ddT, ddL, ddV).  After loading, this creates:
%   x_drv_T_i     = T‐displacement (i from filename)
%   x_drv_L_i     = L‐displacement (i from filename)
%   x_drv_V_i     = V‐displacement (i from filename)
%   ddx_drv_T_i   = T‐acceleration    (i from filename)
%   ddx_drv_L_i   = L‐acceleration    (i from filename)
%   ddx_drv_V_i   = V‐acceleration    (i from filename)
%
% If any expected channel is missing, an error is raised.

  % Split filename into name and extension
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s".', ext);
  end

  % Extract suffix (e.g. 'tgt' or 'DRV') and baseName
  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    baseName = name;
    suffix   = '';
  else
    baseName = name(1:dotIdx-1);
    suffix   = name(dotIdx+1:end);
  end
  suffix = lower(suffix);  % unify to lowercase for comparison

  % Open file
  fid = fopen(filename, 'r');
  if fid < 0
    error('Cannot open "%s" for reading.', filename);
  end

  % Read header lines until 'Data:' marker
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
      % Strip off "Name:" and split on tabs/spaces
      rawNames = strtrim(tline(6:end));
      colNames = regexp(rawNames, '\s+', 'split');
    elseif startsWith(tline, 'No. of Samples:', 'IgnoreCase', true)
      numsamps = sscanf(tline, 'No. of Samples:%f');
    elseif startsWith(tline, 'Time step [s]:', 'IgnoreCase', true)
      dt = sscanf(tline, 'Time step [s]:%f');
    elseif strcmpi(strtrim(tline), 'Data:')
      % Skip the repeated header row (column titles)
      fgetl(fid);
      break;
    end
  end

  if isempty(numsamps) || isempty(dt) || isempty(colNames)
    fclose(fid);
    error('Header did not contain all required fields (Name, No. of Samples, Time step).');
  end

  % Read the numeric block (nCols = numel(colNames), nRows = numsamps)
  ncol = numel(colNames);
  % Build format string like '%f%f%f...' ncol times
  fmt  = repmat('%f', 1, ncol);
  data = fscanf(fid, fmt, [ncol, Inf])';
  fclose(fid);

  if size(data,1) < numsamps
    warning('Read fewer rows (%d) than expected (%d).', size(data,1), numsamps);
    numsamps = size(data,1);
  end

  % Discard any columns that are entirely zero
  nzColsMask = any(data ~= 0, 1);
  data        = data(:, nzColsMask);
  colNames    = colNames(nzColsMask);

  %--- CASE 1: “.tgt” file -----------------------------------------------
  if strcmp(suffix, 'tgt')
    % Expect exactly three non‐zero columns: T, L, V displacements
    % Identify which column corresponds to T, L, V by searching names
    idxT = find(contains(colNames, {'T','t'}, 'IgnoreCase', true), 1);
    idxL = find(contains(colNames, {'L'},     'IgnoreCase', true), 1);
    idxV = find(contains(colNames, {'V'},     'IgnoreCase', true), 1);

    if isempty(idxT) || isempty(idxL) || isempty(idxV)
      error('Could not find T, L, or V columns in .tgt file header.');
    end

    % Compute time vector
    time_vector = (0:(numsamps-1))' * dt;

    % Assign into base workspace
    assignin('base', 'time_vector',   time_vector);
    assignin('base', 'x_tgt_T',       data(:, idxT));
    assignin('base', 'x_tgt_L',       data(:, idxL));
    assignin('base', 'x_tgt_V',       data(:, idxV));

    fprintf('Loaded variables into base workspace:\n');
    fprintf('  time_vector\n');
    fprintf('  x_tgt_T\n');
    fprintf('  x_tgt_L\n');
    fprintf('  x_tgt_V\n');

    return;
  end

  %--- CASE 2: “.drv” file -----------------------------------------------
  if strcmp(suffix, 'drv')
    % Extract numeric index “i” from baseName (assumes format “…_<i>”)
    tokens = regexp(baseName, '_(\d+)$', 'tokens');
    if isempty(tokens)
      error('Could not parse numeric index from baseName "%s".', baseName);
    end
    iStr = tokens{1}{1};  % e.g. '34'

    % Identify columns for displacement (X_T, X_L, X_V) and
    % acceleration (ddX_T, ddX_L, ddX_V).  We do this by matching
    % the header names (case insensitive).
    idxX_T  = find(contains(colNames, {'X_T','x_t','X T'}, 'IgnoreCase', true), 1);
    idxX_L  = find(contains(colNames, {'X_L','x_l','X L'}, 'IgnoreCase', true), 1);
    idxX_V  = find(contains(colNames, {'X_V','x_v','X V'}, 'IgnoreCase', true), 1);
    idxddT  = find(contains(colNames, {'ddX_T','ddx_t','dd X T'}, 'IgnoreCase', true), 1);
    idxddL  = find(contains(colNames, {'ddX_L','ddx_l','dd X L'}, 'IgnoreCase', true), 1);
    idxddV  = find(contains(colNames, {'ddX_V','ddx_v','dd X V'}, 'IgnoreCase', true), 1);

    % Verify all six channels were found
    if any(cellfun(@isempty, {idxX_T, idxX_L, idxX_V, idxddT, idxddL, idxddV}))
      error(['Could not find all required displacement/acceleration columns ', ...
             'in .drv file header. Expected X_T, X_L, X_V, ddX_T, ddX_L, ddX_V.']);
    end

    % Build variable names with index suffix “_iStr”
    var_x_T   = sprintf('x_drv_T_%s',  iStr);
    var_x_L   = sprintf('x_drv_L_%s',  iStr);
    var_x_V   = sprintf('x_drv_V_%s',  iStr);
    var_ddx_T = sprintf('ddx_drv_T_%s',iStr);
    var_ddx_L = sprintf('ddx_drv_L_%s',iStr);
    var_ddx_V = sprintf('ddx_drv_V_%s',iStr);

    % Assign into base workspace
    assignin('base', var_x_T,   data(:, idxX_T));
    assignin('base', var_x_L,   data(:, idxX_L));
    assignin('base', var_x_V,   data(:, idxX_V));
    assignin('base', var_ddx_T, data(:, idxddT));
    assignin('base', var_ddx_L, data(:, idxddL));
    assignin('base', var_ddx_V, data(:, idxddV));

    fprintf('Loaded variables into base workspace:\n');
    fprintf('  %s\n', var_x_T);
    fprintf('  %s\n', var_x_L);
    fprintf('  %s\n', var_x_V);
    fprintf('  %s\n', var_ddx_T);
    fprintf('  %s\n', var_ddx_L);
    fprintf('  %s\n', var_ddx_V);

    return;
  end

  %--- If neither .tgt nor .drv, error out
  error('Unsupported file suffix "%s". Expecting .tgt or .drv (prior to .txt).', suffix);
end

