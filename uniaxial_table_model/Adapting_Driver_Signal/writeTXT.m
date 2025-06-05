function writeTXT(timeVec, dispVec, accVec, filename)
% writeData writes time, displacement, and acceleration data to a text file
%   timeVec  - n-by-1 vector of time values (double)
%   dispVec  - n-by-1 vector of displacement values (double) (PosT)
%   accVec   - n-by-1 vector of acceleration values (double) (accT)
%   filename - string specifying output text file name (e.g., 'output.txt')
%
% The output file will have columns:
% time    PosT    PosL    PosV    accT    accL    accV
% where PosL, PosV, accL, accV are filled with zeros.

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

    % Create zero columns for PosL, PosV, accL, accV
    zeroCol = zeros(nT,1);

    % Combine into one matrix: [time, PosT, PosL, PosV, accT, accL, accV]
    dataMat = [timeVec, dispVec, zeroCol, zeroCol, accVec, zeroCol, zeroCol];

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file %s for writing', filename);
    end

    % Write header line
    fprintf(fid, 'time    PosT    PosL    PosV    accT    accL    accV\n');

    % Write data rows in fixed-point (no scientific notation)
    % Using %.15f prints 15 digits after the decimal point.
    fmt = '%.15f    %.15f    %.15f    %.15f    %.15f    %.15f    %.15f\n';
    for i = 1:nT
        fprintf(fid, fmt, ...
            dataMat(i,1), ... % time
            dataMat(i,2), ... % PosT
            dataMat(i,3), ... % PosL (zero)
            dataMat(i,4), ... % PosV (zero)
            dataMat(i,5), ... % accT
            dataMat(i,6), ... % accL (zero)
            dataMat(i,7));    % accV (zero)
    end

    fclose(fid);
end
