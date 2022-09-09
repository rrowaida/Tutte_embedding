%
% Write a 3D model to a file
%
% function model_write(model, filename);
%
% Input -
%   - model: 3D model structure
%   - filename: name of file to create (the file format is determined
%   from the filename extension)
%
% Output -
%   - status: this variable is 0 if the file was successfully written,
%   or 1 otherwise
%
% See also model
%
function status = model_write(model, filename);
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>
%

    % Select the proper function to write the file, according to the
    % filename extension
    if strcmpi(filename(length(filename)-2:end), 'obj')
        status = model_write_obj(model, filename);
    elseif strcmpi(filename(length(filename)-2:end), 'off')
        status = model_write_off(model, filename);
    else
        status = 1;
        disp(['Could not determine format of file "' filename '" (unknown extension)']);
    end
end
