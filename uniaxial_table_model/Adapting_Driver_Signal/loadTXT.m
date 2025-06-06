function loadTXT(filename)
% loadTXT  Load an LNEC .tgt or .drv text file into the MATLAB workspace,
%           assigning specific variable names based on file type.
%
%   Usage:
%     loadTXT('someFile.tgt.txt')
%       → creates in base workspace:
%            time_vector           (computed from number of samples and dt)
%            x_tgt_T
%            x_tgt_L
%            x_tgt_V
%
%     loadTXT('someFile.drv.txt')
%       → creates in base workspace:
%            x_drv_T
%            x_drv_L
%            x_drv_V
%            ddx_drv_T
%            ddx_drv_L
%            ddx_drv_V
%
%   (If the folder does not contain “.tgt” or “.drv” right before “.txt”,
%   this function will throw an error.)

  %— Split filename into parts (ignore path)
  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s".', ext);
  end

  %— Determine whether this is a “.tgt” or “.drv” file
  %  (look for the substring after the last ‘.’ in “name”)
  dotIdx = find(name == '.', 1, 'last');
  if isempty(dotIdx)
    error('Filename "%s" must contain either “.tgt” or “.drv” before “.txt”.', name);
  end

  baseName = name(1:dotIdx-1);
  suffix   = lower(name(dotIdx+1:end));  % should be either 'tgt' or 'drv'
  if ~ismember(suffix, {'tgt','drv'})
    error('Unrecognized suffix "%s". Expected “tgt” or “drv”.', suffix);
  end

  %— Open the file and read header lines until “Data:”  
  fid = fopen(filename, 'r');
  if fid < 0
    error('Cannot open file "%s".', filename);
  end

  numsamps = [];
  dt       = [];
  while true
    tline = fgetl(fid);
    if ~ischar(tline)
      fclose(fid);
      error('Unexpected end of file while reading header.');
    end

    if startsWith(tline, 'No. of Samples:')
      numsamps = sscanf(tline, 'No. of Samples:%f');
    elseif startsWith(tline, 'Time step [s]:')
      dt = sscanf(tline, 'Time step [s]:%f');
    elseif strcmp(tline, 'Data:')
      fgetl(fid);  % skip the repeated “Name” row
      break;
    end
  end

  if isempty(numsamps) || isempty(dt)
    fclose(fid);
    error('Header did not contain both “No. of Samples” and “Time step [s]”.');
  end

  %— Read the numeric block into an N×M array (rows = samples)
  %   We do not need the column “names” themselves, since we will
  %   assign by position.  However, we need to know how many columns
  %   the file actually has, so we parse the repeated “Name” row.
  %
  %   (At this point, the file pointer is positioned on the first line
  %    of data values.)
  data = fscanf(fid, '%f', [Inf])';
  fclose(fid);

  %— Reshape “data” into numsamps rows by however many columns exist
  %   We know there should be exactly “numsamps” rows, so
  %   number_of_columns = numel(data) / numsamps.
  nTotal       = numel(data);
  nColumns     = nTotal / numsamps;
  if mod(nColumns,1) ~= 0
    error('Data block size (%d) is not consistent with %d samples.', nTotal, numsamps);
  end
  nColumns = floor(nColumns);
  data = reshape(data, nColumns, numsamps)';
  %— Discard any columns that are all zeros
  nzCols = any(data ~= 0, 1);
  data   = data(:, nzCols);
  nCols  = size(data, 2);

  %— Depending on “suffix”, assign to the correct variable names:
  switch suffix
    case 'tgt'
      % Expect exactly 3 nonzero columns: [x_tgt_T, x_tgt_L, x_tgt_V].
      if nCols < 3
        error('“.tgt” file must contain at least 3 nonzero columns. Found %d.', nCols);
      end

      % Assign displacements:
      x_tgt_T = data(:, 1);
      x_tgt_L = data(:, 2);
      x_tgt_V = data(:, 3);

      % Build time_vector (numsamps × 1):
      time_vector = (0:numsamps-1)' * dt;

      % Push into base workspace:
      assignin('base', 'x_tgt_T',    x_tgt_T);
      assignin('base', 'x_tgt_L',    x_tgt_L);
      assignin('base', 'x_tgt_V',    x_tgt_V);
      assignin('base', 'time_vector', time_vector);

      loadedVars = { 'time_vector', 'x_tgt_T', 'x_tgt_L', 'x_tgt_V' };

    case 'drv'
      % Expect exactly 6 nonzero columns:
      %   [x_drv_T, x_drv_L, x_drv_V, ddx_drv_T, ddx_drv_L, ddx_drv_V].
      if nCols < 6
        error('“.drv” file must contain at least 6 nonzero columns. Found %d.', nCols);
      end

      % Assign first three to displacement, next three to acceleration:
      x_drv_T   = data(:, 1);
      x_drv_L   = data(:, 2);
      x_drv_V   = data(:, 3);
      ddx_drv_T = data(:, 4);
      ddx_drv_L = data(:, 5);
      ddx_drv_V = data(:, 6);

      % Push into base workspace:
      assignin('base', 'x_drv_T',    x_drv_T);
      assignin('base', 'x_drv_L',    x_drv_L);
      assignin('base', 'x_drv_V',    x_drv_V);
      assignin('base', 'ddx_drv_T',  ddx_drv_T);
      assignin('base', 'ddx_drv_L',  ddx_drv_L);
      assignin('base', 'ddx_drv_V',  ddx_drv_V);

      loadedVars = { 'x_drv_T', 'x_drv_L', 'x_drv_V', ...
                     'ddx_drv_T', 'ddx_drv_L', 'ddx_drv_V' };

    otherwise
      % (This branch should never be reached, since we checked “suffix” above.)
      error('Internal error: unrecognized suffix "%s".', suffix);
  end

  %— Print out what we loaded
  fprintf('Loaded variables into base workspace:\n');
  for k = 1:numel(loadedVars)
    fprintf('  %s\n', loadedVars{k});
  end
end

