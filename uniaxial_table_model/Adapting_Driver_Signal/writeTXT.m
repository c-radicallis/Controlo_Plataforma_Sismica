function writeTXT(timeVec, dispVec, accVec, filename)
% writeTXT writes time, displacement, and acceleration data to a text file
%   timeVec  - n-by-1 vector of time values (numeric, will be cast to int32)
%   dispVec  - n-by-1 vector of displacement values (numeric, will be cast to int32)
%   accVec   - n-by-1 vector of acceleration values (numeric, will be cast to int32)
%   filename - string specifying output text file name (e.g., 'output.txt')
%
% The output file will have columns (all integer, no scientific notation):
%   time   PosT   PosL   PosV   accT   accL   accV
% where PosL, PosV, accL, accV are filled with zeros (int32).

    % Ensure inputs are column vectors of same length
    nT = numel(timeVec);
    nD = numel(dispVec);
    nA = numel(accVec);
    if nT ~= nD || nT ~= nA
        error('timeVec, dispVec, and accVec must have the same length');
    end
    timeVec = timeVec(:);
    dispVec = dispVec(:);
    accVec  = accVec(:);

    % Cast to int32 (LabVIEW signed 32-bit)
    timeInt = int32(timeVec);
    dispInt = int32(dispVec);
    accInt  = int32(accVec);

    % Create zero columns (int32)
    zeroCol = int32(zeros(nT,1));

    % Combine into one matrix of int32: [time, PosT, PosL, PosV, accT, accL, accV]
    dataMat = [ timeInt, dispInt, zeroCol, zeroCol, accInt, zeroCol, zeroCol ];

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file %s for writing', filename);
    end

    % Write header line (note: header is text, LabVIEW will skip or parse it)
    fprintf(fid, 'time  PosT  PosL  PosV  accT  accL  accV\n');

    % Integer format: %d for each of the 7 columns, separated by spaces
    fmt = '%d  %d  %d  %d  %d  %d  %d\n';
    for i = 1:nT
        fprintf(fid, fmt, ...
            dataMat(i,1), dataMat(i,2), dataMat(i,3), ...
            dataMat(i,4), dataMat(i,5), dataMat(i,6), ...
            dataMat(i,7) );
    end

    fclose(fid);
end
