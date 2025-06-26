function writeTXT_then_LTF(t_vector , x , ddx , folder , filename)
%  "Folder" is the directory from which the file is read and writen

    writeTXT(t_vector , x , ddx , folder , filename)

    full_filename = fullfile(folder, filename);
    
    try
        py_output = py.TXT_to_LTF.txt_to_ltf(full_filename,folder);
        % py_output is a Python string; convert to MATLAB char:
        output_path = char(py_output);
        fprintf('Python function returned output path: %s\n', output_path);
    catch ME
        disp('Error calling Python function:');
        disp(ME.message);
    end

end