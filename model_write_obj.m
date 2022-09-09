%
% Write a 3D model to an obj file
%
% function model_write_obj(model, filename);
%
% Input -
%   - model: 3D model structure
%   - filename: name of obj file to create
%
% Output -
%   - status: this variable is 0 if the file was successfully written,
%   or 1 otherwise
%
% The function writes the vertices and faces of the model to an obj
% file. If vertex colors are present, these are also written into the
% file following the extended obj format.
%
% See also model, model_write
%
function status = model_write_obj(model, filename);
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>
%

    % Pre-process colors, if defined
    vertex_color = 0;
    face_color = 0; % Not used now, but can be used in the future
    clr = [];
    if (isfield(model.mesh, 'FaceVertexCData'))
        % Check type of color
        if size(model.mesh.FaceVertexCData, 1) == size(model.mesh.vertices, 1)
            vertex_color = 1;
        elseif size(model.mesh.FaceVertexCData, 1) == size(model.mesh.face, 1)
            face_color = 1;
        end
        % Check if colors are defined as RGB or a function
        if vertex_color == 1%|| face_color
            if size(model.mesh.FaceVertexCData, 2) == 1
                % Map function values to RGB data
                % Map function values from [min, max] to [0, 2/3], which is
                % red to blue hue in an HSV color map
                func = map_val(model.mesh.FaceVertexCData, 0, 2/3);
                % Use max(func)-func as the hue value and get RGB colors
                % Saturation and brightness are set to 0.8
                clr = hsv2rgb([max(func)-func ones(length(func), 1)*0.8 ...
                                    ones(length(func),1)*0.8]);
            else
                % If we have RGB data, just use it
                clr = model.mesh.FaceVertexCData;
            end
        end
    end

    % Open file
    fid = fopen(filename, 'w');
    status = 0;
    if fid == -1
        disp(['ERROR: could not open file "' filename '"']);
        status = 1;
        return;
    end

    % Write vertices
    if vertex_color
        for i = 1:size(model.mesh.vertices, 1)
            fprintf(fid, 'v %f %f %f %f %f %f\n', ...
                        model.mesh.vertices(i, 1), model.mesh.vertices(i, 2), model.mesh.vertices(i, 3), ...
                        clr(i, 1), clr(i, 2), clr(i, 3));
        end
    else
        for i = 1:size(model.mesh.vertices, 1)
            fprintf(fid, 'v %f %f %f\n', ...
                        model.mesh.vertices(i, 1), model.mesh.vertices(i, 2), model.mesh.vertices(i, 3));

        end
    end

    % Write faces
    for i = 1:size(model.mesh.faces, 1)
        fprintf(fid, 'f %d %d %d\n', ...
                    model.mesh.faces(i, 1), model.mesh.faces(i, 2), model.mesh.faces(i, 3));
    end

    % Close file
    fclose(fid);
end
