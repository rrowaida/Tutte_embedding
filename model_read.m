%
% Read a 3D model from a file
%
% function [model, status] = model_read(filename)
%
% Input -
%   - filename: name of the file to load (the file format is determined
%   from the filename extension)
%
% Output -
%   - model: 3D model structure
%   - status: this variable is 0 if the file was successfully opened, or
%   1 otherwise
%
% See also model
%
function [model, status] = model_read(filename)
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>
%

    % Select the proper function to read the file, according to the filename
    % extension
    if strcmpi(filename(length(filename)-2:end), 'obj')
        [model, status] = model_read_obj(filename);
    elseif strcmpi(filename(length(filename)-2:end), 'off')
        [model, status] = model_read_off(filename);
    else
        model = [];
        status = 1;
        disp(['Could not determine format of file "' filename '" (unknown extension)']);
    end
end
