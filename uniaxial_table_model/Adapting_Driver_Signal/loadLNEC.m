function loadLNEC(filename)
% LOADLNEC  Load an LNEC .tgt or .drv text file into MATLAB workspace,
% prefixing each variable with the file’s suffix (e.g. 'tgtDispT').
%
%   loadLNEC('LAquilaReducedScale.tgt.txt')
%   ⇒ creates variables tgtDispT, tgtDispL, ..., tgtt

  [~, name, ext] = fileparts(filename);
  if ~strcmpi(ext, '.txt')
    error('Expected a .txt file. Got "%s"', ext);
  end

  % extract suffix after last '.' in the name
  dotIdx = find(name=='.', 1, 'last');
  if isempty(dotIdx)
    prefix = matlab.lang.makeValidName(name);
  else
    prefix = matlab.lang.makeValidName(name(dotIdx+1:end));
  end

  %--- Open file
  fid = fopen(filename,'r');
  if fid<0, error('Cannot open %s', filename); end

  %--- Read header
  while true
    tline = fgetl(fid);
    if ~ischar(tline), error('Unexpected EOF in header.'); end
    if startsWith(tline,'Name:')
      names = strsplit(strtrim(tline(6:end)), '\t');
    elseif startsWith(tline,'No. of Samples:')
      numsamps = sscanf(tline, 'No. of Samples:%f');
    elseif startsWith(tline,'Time step [s]:')
      dt = sscanf(tline, 'Time step [s]:%f');
    elseif strcmp(tline,'Data:')
      % skip the repeated header line
      fgetl(fid);
      break
    end
  end

  %--- Read the data
  ncol = numel(names);
  data = fscanf(fid, repmat('%f',1,ncol), [ncol Inf])';
  fclose(fid);

  %--- Discard all‑zero columns
  nz = any(data~=0,1);
  data = data(:,nz);
  names = names(nz);

  %--- Assign to workspace with prefix
  for i = 1:numel(names)
    var = matlab.lang.makeValidName(names{i});
    fullname = [prefix var];
    assignin('base', fullname, data(:,i));
  end

  %--- Time vector
  t = (0:numsamps-1)' * dt;
  assignin('base','time_vector', t);

  fprintf('Loaded %d channels into workspace with prefix "%s"; created %s (time vector).\n', ...
          numel(names), prefix, [prefix 't']);
end

