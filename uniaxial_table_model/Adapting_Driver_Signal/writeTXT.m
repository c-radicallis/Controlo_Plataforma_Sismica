function writeTXT(timeVec, dispVec, accVec, folderPath, filename)
% writeTXT writes time, displacement, and acceleration data to a text file
%   timeVec    - n×1 vector of time values (double)
%   dispVec    - n×1 vector of displacement values (double) (PosT)
%   accVec     - n×1 vector of acceleration values (double) (accT)
%   folderPath - string specifying the folder in which to save the file
%   filename   - string specifying output text file name (e.g., 'output' or 'output.txt')
%
% The output file will have columns:
% time    PosT    PosL    PosV    accT    accL    accV
% where PosL, PosV, accL, accV are filled with zeros.
%
% If filename does not already end in “.txt” (case-insensitive), “.txt” is appended.    

    % If the folder does not exist, attempt to create it
    if ~exist(folderPath, 'dir')
        mkdirStatus = mkdir(folderPath);
        if ~mkdirStatus
            error('Could not create directory: %s', folderPath);
        end
    end

    % Ensure filename ends with '.txt' (case-insensitive)
    if ~endsWith(filename, '.txt', 'IgnoreCase', true)
        filename = [filename, '.txt'];
    end

    % Build the full filename (folder + filename)
    fullPath = fullfile(folderPath, filename);

    % Combine into one matrix: [time, PosT, PosL, PosV, accT, accL, accV]
    dim = size(dispVec);
    dim1= dim(1);
    % Create zero columns 
    zeroCol = zeros(dim1,1);

    if dim(2)==2
        dataMat = [timeVec, dispVec(:,1) , dispVec(:,2) , zeroCol, accVec(:,1), accVec(:,2), zeroCol];
    else
        dataMat = [timeVec, dispVec , zeroCol , zeroCol, accVec , zeroCol , zeroCol];
    end

    % Open file for writing
    fid = fopen(fullPath, 'w');
    if fid == -1
        error('Could not open file %s for writing', fullPath);
    end

    % Write header line
    fprintf(fid, 'time    PosT    PosL    PosV    accT    accL    accV\n');

    % Write data rows in fixed-point (no scientific notation)
    % Using %.3f for time and %.18f for others, adjust as needed
    fmt = '%.3f    %.18f    %.18f    %.18f    %.18f    %.18f    %.18f\n';
    for i = 1:dim1
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
