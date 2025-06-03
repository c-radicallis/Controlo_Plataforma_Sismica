function writeTXT(timeVec, dispVec, accVec, filename)
% writeData writes time, displacement, and acceleration data to a text file
%   timeVec  - n-by-1 vector of time values (double)
%   dispVec  - n-by-1 vector of displacement values (double) (PosT)
%   accVec   - n-by-1 vector of acceleration values (double) (accT)
%   filename - string specifying output text file name (e.g., 'output.txt')
%
% The output file will have columns:
% time  PosT  PosL  PosV  accT  accL  accV
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
    accVec = accVec(:);

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
    fprintf(fid, 'time  PosT  PosL  PosV  accT  accL  accV\n');

    % Write data rows with full double precision
    % Using %.15g to preserve up to 15 significant digits
    fmt = '%.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g\n';
    for i = 1:nT
        fprintf(fid, fmt, dataMat(i,1), dataMat(i,2), dataMat(i,3), dataMat(i,4), dataMat(i,5), dataMat(i,6), dataMat(i,7));
    end

    fclose(fid);
end