function loadTXT(filename)
% loadTXT  Load an LNEC .tgt, .drv or .acq text file into MATLAB workspace,
%          creating appropriately named variables in the base workspace.
%          Ignores columns of zeros or containing NaNs.

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
      colNames = strsplit(rawNames, '\t');
    elseif startsWith(tline, 'Type:', 'IgnoreCase', true)
      rawTypes = strtrim(tline(6:end));
      typeNames = strsplit(rawTypes, '\t');
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

  if isempty(numsamps) || isempty(dt) || isempty(colNames)
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
  % 4b. Discard columns all-zero or containing NaN %
  %--------------------------------------------%
  validCols = any(data ~= 0, 1) & ~any(isnan(data), 1);
  data      = data(:, validCols);
  colNames  = colNames(validCols);
  if ~isempty(typeNames)
    typeNames = typeNames(validCols);
  end

  %--------------------------------------------%
  % 5. Decide behavior based on file suffix    %
  %--------------------------------------------%
  switch suffix
    case 'tgt'
      % ensure types trimmed
      typeNames  = typeNames;
      time_vector = (0:(numsamps-1))' * dt;
      assignin('base', 'time_vector', time_vector);
      fprintf('Loaded: time_vector\n');

      dispIdxs = find(strcmpi(typeNames, 'Displacement'));
      accIdxs  = find(strcmpi(typeNames, 'Acceleration'));

      idxDispT = findChannelIndex(colNames, dispIdxs, 'T');
      idxDispL = findChannelIndex(colNames, dispIdxs, 'L');
      idxDispV = findChannelIndex(colNames, dispIdxs, 'V');
      idxAccT  = findChannelIndex(colNames, accIdxs,  'T');
      idxAccL  = findChannelIndex(colNames, accIdxs,  'L');
      idxAccV  = findChannelIndex(colNames, accIdxs,  'V');

      if ~isempty(idxDispT), assignin('base','x_tgt_T',data(:,idxDispT));fprintf('Loaded: x_tgt_T\n'); end
      if ~isempty(idxDispL), assignin('base','x_tgt_L',data(:,idxDispL));fprintf('Loaded: x_tgt_L\n'); end
      if ~isempty(idxDispV), assignin('base','x_tgt_V',data(:,idxDispV));fprintf('Loaded: x_tgt_V\n'); end
      if ~isempty(idxAccT),  assignin('base','ddx_tgt_T',data(:,idxAccT));fprintf('Loaded: ddx_tgt_T\n'); end
      if ~isempty(idxAccL),  assignin('base','ddx_tgt_L',data(:,idxAccL));fprintf('Loaded: ddx_tgt_L\n'); end
      if ~isempty(idxAccV),  assignin('base','ddx_tgt_V',data(:,idxAccV));fprintf('Loaded: ddx_tgt_V\n'); end
      return

    case 'drv'
      tokens = regexp(baseName, '_(\d+)$', 'tokens');
      if isempty(tokens)
        error('Cannot parse index from filename "%s".', baseName);
      end
      iStr = tokens{1}{1};
      time_drv = (0:(numsamps-1))' * dt;
      assignin('base', sprintf('time_drv_%s',iStr), time_drv);

      dispIdxs = find(strcmpi(typeNames, 'Displacement'));
      if isempty(dispIdxs), dispIdxs = 1:numel(typeNames); end

      idxT = findChannelIndex(colNames, dispIdxs, 'T');
      idxL = findChannelIndex(colNames, dispIdxs, 'L');
      idxV = findChannelIndex(colNames, dispIdxs, 'V');

      if ~isempty(idxT), assignin('base',sprintf('x_drv_T_%s',iStr),data(:,idxT)); end
      if ~isempty(idxL), assignin('base',sprintf('x_drv_L_%s',iStr),data(:,idxL)); end
      if ~isempty(idxV), assignin('base',sprintf('x_drv_V_%s',iStr),data(:,idxV)); end
      return

    case 'acq'
      % Treat like .tgt but name variables with _acq suffix
      time_acq = (0:(numsamps-1))' * dt;
      assignin('base', 'time_acq', time_acq);
      fprintf('Loaded: time_acq\n');

      dispIdxs = find(strcmpi(typeNames, 'Displacement'));
      accIdxs  = find(strcmpi(typeNames, 'Acceleration'));

      idxDispT = findChannelIndex(colNames, dispIdxs, 'T');
      idxDispL = findChannelIndex(colNames, dispIdxs, 'L');
      idxDispV = findChannelIndex(colNames, dispIdxs, 'V');
      idxAccT  = findChannelIndex(colNames, accIdxs, 'T');
      idxAccL  = findChannelIndex(colNames, accIdxs, 'L');
      idxAccV  = findChannelIndex(colNames, accIdxs, 'V');

      if ~isempty(idxDispT), assignin('base','x_acq_T',data(:,idxDispT));fprintf('Loaded: x_acq_T\n'); end
      if ~isempty(idxDispL), assignin('base','x_acq_L',data(:,idxDispL));fprintf('Loaded: x_acq_L\n'); end
      if ~isempty(idxDispV), assignin('base','x_acq_V',data(:,idxDispV));fprintf('Loaded: x_acq_V\n'); end
      if ~isempty(idxAccT),  assignin('base','ddx_acq_T',data(:,idxAccT));fprintf('Loaded: ddx_acq_T\n'); end
      if ~isempty(idxAccL),  assignin('base','ddx_acq_L',data(:,idxAccL));fprintf('Loaded: ddx_acq_L\n'); end
      if ~isempty(idxAccV),  assignin('base','ddx_acq_V',data(:,idxAccV));fprintf('Loaded: ddx_acq_V\n'); end
      return

    otherwise
      error('Unsupported suffix "%s". Expect tgt, drv, or acq.', suffix);
  end
end

function idx = findChannelIndex(colNames, idxList, letter)
  idx = [];
  for jj = idxList
    if contains(colNames{jj}, letter, 'IgnoreCase', true)
      idx = jj;
      return;
    end
  end
end
