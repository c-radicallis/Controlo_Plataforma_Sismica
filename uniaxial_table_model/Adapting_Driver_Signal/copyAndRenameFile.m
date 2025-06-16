function status = copyAndRenameFile(srcPath, destDir, newName)
%COPYANDRENAMEFILE Copies a file to a new location and renames it.
%   status = COPYANDRENAMEFILE(srcPath, destDir, newName) attempts to copy
%   the file at full path srcPath into the directory destDir with filename
%   newName. If destDir does not exist, it will be created. 
%
%   Inputs:
%     srcPath  - Full path to the existing source file, e.g.:
%                'C:\folder\oldname.txt' or '/home/user/oldname.txt'
%     destDir  - Path to the folder where you want to copy, e.g.:
%                'D:\backup' or '/home/user/backup'
%     newName  - New filename (with extension) in the destination, e.g.:
%                'newname.txt'
%
%   Output:
%     status   - Logical true if copy succeeded, false otherwise.
%
%   Example:
%     ok = copyAndRenameFile('C:\data\report.pdf', 'D:\archive', 'report_old.pdf');
%     if ok
%         fprintf('Copy succeeded.\n');
%     else
%         fprintf('Copy failed.\n');
%     end

    % Initialize status
    status = false;

    % Validate srcPath
    if ~ischar(srcPath) && ~isstring(srcPath)
        error('srcPath must be a character vector or string.');
    end
    srcPath = char(srcPath);

    if exist(srcPath, 'file') ~= 2
        error('Source file does not exist: %s', srcPath);
    end

    % Validate destDir
    if ~ischar(destDir) && ~isstring(destDir)
        error('destDir must be a character vector or string.');
    end
    destDir = char(destDir);

    % If destDir does not exist, try to create it
    if exist(destDir, 'dir') ~= 7
        try
            mkdir(destDir);
            fprintf('Destination directory did not exist; created: %s\n', destDir);
        catch ME
            error('Failed to create destination directory "%s": %s', destDir, ME.message);
        end
    end

    % Validate newName
    if ~ischar(newName) && ~isstring(newName)
        error('newName must be a character vector or string.');
    end
    newName = char(newName);
    % Optionally, you could check that newName has no path separators:
    [~, nameOnly, ext] = fileparts(newName);
    if isempty(nameOnly)
        error('newName must contain a filename (with or without extension).');
    end
    % If extension missing, you may choose to warn or append same as source:
    % E.g., if isempty(ext), ext = fileparts(srcPath) extension; but here we trust user.

    % Construct full destination path
    destPath = fullfile(destDir, [nameOnly ext]);

    % If destPath is the same as srcPath (e.g., same folder and same name), warn or error
    fullSrc = which(absPath(srcPath)); %#ok<ASGLU> % normalize path if possible
    % But which() may not resolve arbitrary paths; skip strict check. Instead:
    if strcmpi(fullfile(srcPath), fullfile(destPath))
        warning('Source and destination paths are identical. No action taken.');
        status = false;
        return;
    end

    % Attempt to copy
    try
        [copySuccess, msg, msgID] = copyfile(srcPath, destPath);
        if copySuccess
            fprintf('File copied successfully to:\n    %s\n', destPath);
            status = true;
        else
            warning('copyfile failed: %s (ID: %s)', msg, msgID);
            status = false;
        end
    catch ME
        warning('An error occurred while copying: %s', ME.message);
        status = false;
    end
end

function p = absPath(p)
% ABSPATH Return an absolute path if possible (best-effort).
%   This helper tries to resolve relative paths to full paths.
    if isfolder(fileparts(p))
        % If p includes folder components, make absolute relative to pwd
        if ~isAbsolutePath(p)
            p = fullfile(pwd, p);
        end
    else
        % If file itself may not exist in MATLAB path, just attempt to resolve
        if ~isAbsolutePath(p)
            p = fullfile(pwd, p);
        end
    end
end

function tf = isAbsolutePath(p)
% ISABSOLUTEPATH Simple check for absolute path (Windows or Unix)
    if ispc
        % Windows: absolute if starts with drive letter, e.g. 'C:\' or '\\' for UNC
        tf = ~isempty(regexp(p, '^[A-Za-z]:\\', 'once')) || strncmp(p, '\\', 2);
    else
        % Unix: starts with '/'
        tf = strncmp(p, '/', 1);
    end
end
