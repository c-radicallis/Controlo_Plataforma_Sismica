function writeTXT_then_LTF(t_vector , x_acq , ddx_acq , save_folder , filename_acq)
    
    writeTXT(t_vector , x_acq , ddx_acq , save_folder , filename_acq)
    full_filename_acq = fullfile(save_folder, filename_acq);
    try
        py_output = py.TXT_to_LTF.txt_to_ltf(full_filename_acq,save_folder);
        % py_output is a Python string; convert to MATLAB char:
        output_path = char(py_output);
        fprintf('Python function returned output path: %s\n', output_path);
    catch ME
        disp('Error calling Python function:');
        disp(ME.message);
    end

end