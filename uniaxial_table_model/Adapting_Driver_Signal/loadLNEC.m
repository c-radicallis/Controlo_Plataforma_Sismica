function loadLNEC(filename)
% LOADLNEC  Load an LNEC .tgt or .drv text file into MATLAB workspace,
% prints the names of loaded variables.
%
% Usage:
%   loadLNEC('LAquilaReducedScale.tgt.txt')
%
% Reads the header to determine channel names, number of samples, and time step.
% Loads the data block, discards columns that are entirely zero, then assigns each
% non-null column to the base workspace with a prefix based on the filename suffix.
% Also creates a time vector "time_vector" (length = No. of Samples, spacing = Time step).
% Finally, prints the names of all variables it loaded into the workspace.

  % Split filename into parts
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s"', ext);
  end

  % Extract suffix after last '.' in the (basename) to use as prefix
  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    % If no dot in name, use entire name as prefix (made valid)
    prefix = matlab.lang.makeValidName(name);
  else
    % Take the part after last dot (e.g. 'tgt' from 'LAquila.tgt')
    prefix = matlab.lang.makeValidName(name(dotIdx+1:end));
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
    % Capture channel names
    if startsWith(tline, 'Name:')
      % Everything after 'Name:' is a tab-delimited list of channel labels
      names = strsplit(strtrim(tline(6:end)), '\t');
    % Capture number of samples
    elseif startsWith(tline, 'No. of Samples:')
      numsamps = sscanf(tline, 'No. of Samples:%f');
    % Capture time step (in seconds)
    elseif startsWith(tline, 'Time step [s]:')
      dt = sscanf(tline, 'Time step [s]:%f');
    % When we see 'Data:', the next line is the repeated header, so skip it
    elseif strcmp(tline, 'Data:')
      fgetl(fid);  % skip the repeated "Name" row
      break;
    end
  end

  %--- Read the numeric data block into an N-by-M array (rows=samples)
  ncol = numel(names);  % how many columns to expect per row
  data = fscanf(fid, repmat('%f', 1, ncol), [ncol, Inf])';
  fclose(fid);

  %--- Identify and discard columns that are all zeros
  nzCols = any(data ~= 0, 1);      % logical mask of non-zero columns
  data   = data(:, nzCols);        % keep only non-null columns
  names  = names(nzCols);          % matching subset of channel names

  %--- Assign each non-null column into the base workspace
  loadedVars = {};  % track for printing
  for i = 1:numel(names)
    varLabel = matlab.lang.makeValidName(names{i});
    fullName = [prefix varLabel];
    assignin('base', fullName, data(:, i));
    loadedVars{end+1} = fullName;  %#ok<AGROW>
  end

  %--- Build and assign the time vector (common for all channels)
  time_vector = (0:numsamps-1)' * dt;  % column vector of length numsamps
  assignin('base', 'time_vector', time_vector);
  loadedVars{end+1} = 'time_vector';

  %--- Print the names of loaded variables
  fprintf('Loaded variables into base workspace:\n');
  for k = 1:numel(loadedVars)
    fprintf('  %s\n', loadedVars{k});
  end
end
