function loadTXT(filename)
% LOADLNEC  Load an LNEC .tgt or .drv text file into MATLAB workspace,
% printing the names of loaded variables, formatted as suffix+channel+'_baseName'.
%
% Usage:
%   loadLNEC('LAquilaReducedScale.tgt.txt')
%   → creates variables like 'tgtDispT_LAquilaReducedScale' and 'tgt_time_vector_LAquilaReducedScale'.

  % Split filename into parts (path is ignored)
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s"', ext);
  end

  % Extract suffix after last '.' and baseName before that
  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    % If no dot in name, suffix is empty, baseName is entire name
    baseName = matlab.lang.makeValidName(name);
    suffix   = '';
  else
    baseName = matlab.lang.makeValidName(name(1:dotIdx-1));
    suffix   = matlab.lang.makeValidName(name(dotIdx+1:end));
  end

  %--- Open the file for reading
  fid = fopen(filename, 'r');
  if fid < 0
    error('Cannot open %s', filename);
  end

  %--- Read header lines until 'Data:' marker
  while true
    tline = fgetl(fid);
    if ~ischar(tline)
      error('Unexpected end of file while reading header.');
    end
    if startsWith(tline, 'Name:')
      names = strsplit(strtrim(tline(6:end)), '\t');
    elseif startsWith(tline, 'No. of Samples:')
      numsamps = sscanf(tline, 'No. of Samples:%f');
    elseif startsWith(tline, 'Time step [s]:')
      dt = sscanf(tline, 'Time step [s]:%f');
    elseif strcmp(tline, 'Data:')
      fgetl(fid);  % skip the repeated "Name" row
      break;
    end
  end

  %--- Read the numeric data block into an N-by-M array (rows=samples)
  ncol = numel(names);
  data = fscanf(fid, repmat('%f', 1, ncol), [ncol, Inf])';
  fclose(fid);

  %--- Identify and discard columns that are all zeros
  nzCols = any(data ~= 0, 1);
  data   = data(:, nzCols);
  names  = names(nzCols);

  %--- Assign each non-null column into the base workspace
  loadedVars = {};
  for i = 1:numel(names)
    varLabel = matlab.lang.makeValidName(names{i});
    if isempty(suffix)
      fullName = sprintf('%s_%s', varLabel, baseName);
    else
      fullName = sprintf('%s%s_%s', suffix, varLabel, baseName);
    end
    assignin('base', fullName, data(:, i));
    loadedVars{end+1} = fullName;  %#ok<AGROW>
  end

  %--- Build and assign the time vector with the same naming convention
  if isempty(suffix)
    timeVectorName = sprintf('time_vector_%s', baseName);
  else
    timeVectorName = sprintf('%s_time_vector_%s', suffix, baseName);
  end
  time_vector = (0:numsamps-1)' * dt;
  assignin('base', timeVectorName, time_vector);
  loadedVars{end+1} = timeVectorName;

  %--- Print the names of loaded variables
  fprintf('Loaded variables into base workspace:\n');
  for k = 1:numel(loadedVars)
    fprintf('  %s\n', loadedVars{k});
  end
end
