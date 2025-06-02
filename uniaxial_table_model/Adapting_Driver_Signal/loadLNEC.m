function loadLNEC(filename)
% load_tgt.m
% Script to read an LNEC .tgt file, discard null columns, and load
% non-null data as workspace vectors prefixed with the file base name.

% Specify filename
filename = 'LAquilaReducedScale.tgt';
[~, base, ~] = fileparts(filename);

fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open file: %s', filename);
end

% Read header lines until 'Data:' marker
headerLines = {};
while true
    tline = fgetl(fid);
    if ~ischar(tline)
        error('Unexpected end of file before Data section');
    end
    headerLines{end+1} = strtrim(tline);
    if strcmpi(strtrim(tline), 'Data:')
        break;
    end
end

% Next line: column names	names separated by tabs
colNames = strsplit(strtrim(fgetl(fid)), '\t');
numCols = numel(colNames);

% Parse 'No. of Samples:' and 'Time step [s]:' from headerLines
% Find line indices
nSamplesLine = find(startsWith(headerLines, 'No. of Samples:'), 1);
dtLine       = find(startsWith(headerLines, 'Time step [s]:'), 1);

% Extract numeric lists
samplesList = sscanf(headerLines{nSamplesLine}, 'No. of Samples:\s*%f', 1);
dtValue     = sscanf(headerLines{dtLine},       'Time step [s]:\s*%f', 1);

% Read data matrix: numCols columns, unknown rows
fmt = repmat('%f', 1, numCols);
dataArray = textscan(fid, fmt, 'Delimiter', '\t');
fclose(fid);

% Convert cell to matrix: rows x cols
dataMat = cell2mat(dataArray);

% Identify and discard null (all-zero) columns
nullCols = all(dataMat == 0, 1);
keepCols = find(~nullCols);

% Create time vector
N = samplesList;
timeVec = (0:N-1) * dtValue;

% Assign non-null columns to base workspace with names "<base><colName>"
for ii = keepCols
    varname = matlab.lang.makeValidName([base colNames{ii}]);
    assignin('base', varname, dataMat(:, ii));
end

% Assign time vector as <base>Time
timeVar = matlab.lang.makeValidName([base 'Time']);
assignin('base', timeVar, timeVec);

fprintf('Loaded %d columns into workspace, created variables prefixed ''%s'' and time vector ''%s''.\n', ...
        numel(keepCols), base, timeVar);

end

