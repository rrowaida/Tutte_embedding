%
% Read an obj file and return a 3D model structure
%
% function [model, status] = model_read_obj(filename)
%
% Input -
%   - filename: name of obj file to load
%
% Output -
%   - model: 3D model structure
%   - status: this variable is 0 if the file was successfully opened, or
%   1 otherwise
%
% The function reads the vertices and faces of the mesh, and ignores
% other information such as normals and texture coordinates
%
% See also model, model_read
%
function [model, status] = model_read_obj(filename)
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>, Diego Gonzalez
%

    % Open file
    fid = fopen(filename, 'r');
    status = 0;
    if fid == -1
        disp(['ERROR: could not open file "' filename '"']);
        model = [];
        status = 1;
        return;
    end

    % Index of current vertex, face, and color to be added
    vindex = 1;
    findex = 1;
    cindex = 1;

    % Initialize data matrices
    X = zeros(0, 3);
    F = zeros(0, 3);

    % Read content
    while(feof(fid) ~= 1)
        % Read one line
        line = '';
        line = fgetl(fid);
        % A line needs at least two characters to be meaningful
        if length(line) <= 1
            continue;
        end
        % Parse line according to obj format
        % Vertex
        if (line(1) == 'v') && (line(2) == ' ')
            line = line(3:length(line));
            X(vindex, :) = sscanf(line, '%f %f %f');
            vindex = vindex + 1;
        % Face, which can come in many different formats
        elseif (line(1) == 'f') && (line(2) == ' ')
            % Get data portion of line
            line = line(3:length(line));
            % Try to read quad
            processed = 0;
            % Quad: first format
            temp = sscanf(line, '%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');
            if length(temp) >= 12
                verts = temp([1 4 7 10]);
                F(findex, :) = verts([1 2 3]);
                findex = findex + 1;
                F(findex, :) = verts([1 3 4]);
                processed = 1;
            end
            if ~processed
                % Quad: second format
                temp = sscanf(line, '%d//%d %d//%d %d//%d %d//%d');
                if length(temp) >= 8
                    % If it worked, create two triangles from the quad
                    verts = temp([1 3 5 7]);
                    F(findex, :) = verts([1 2 3]);
                    findex = findex + 1;
                    F(findex, :) = verts([1 3 4]);
                    processed = 1;
                end
            end
            if ~processed
                % Quad: third format
                temp = sscanf(line, '%d/%d %d/%d %d/%d %d/%d');
                if length(temp) >= 8
                    % If it worked, create two triangles from the quad
                    verts = temp([1 3 5 7]);
                    F(findex, :) = verts([1 2 3]);
                    findex = findex + 1;
                    F(findex, :) = verts([1 3 4]);
                    processed = 1;
                end
            end
            if ~processed
                % Quad: fourth format
                temp = sscanf(line, '%d %d %d %d');
                if length(temp) >= 4
                    % If it worked, create two triangles from the quad
                    verts = temp([1 2 3 4]);
                    F(findex, :) = verts([1 2 3]);
                    findex = findex + 1;
                    F(findex, :) = verts([1 3 4]);
                    processed = 1;
                end
            end
            % Try to read triangle
            if ~processed
                % Triangle: first format
                temp = sscanf(line, '%d/%d/%d %d/%d/%d %d/%d/%d');
                if length(temp) >= 9
                    F(findex, :) = temp([1 4 7]);
                    processed = 1;
                end
            end
            if ~processed
                % Triangle: second format
                temp = sscanf(line, '%d//%d %d//%d %d//%d');
                if length(temp) >= 6
                    % Discard the extra info and only use vertex refs
                    F(findex, :) = temp([1 3 5]);
                    processed = 1;
                end
            end
            if ~processed
                % Triangle: third format
                temp = sscanf(line, '%d/%d %d/%d %d/%d');
                if length(temp) >= 6
                    % Discard the extra info and only use vertex refs
                    F(findex, :) = temp([1 3 5]);
                    processed = 1;
                end
            end
            if ~processed
                % Triangle: fourth format
                F(findex, :) = sscanf(line, '%d %d %d');
            end
            findex = findex + 1;
        % Color
        elseif strcmp(line(1:2), 'vc') || strcmp(line(1:2), 'fc')
            line = line(3:length(line));
            C(cindex, :) = sscanf(line, '%f %f %f');
            cindex = cindex + 1;
        end
    end

    % Close file
    fclose(fid);

    % Set up output model
    model = struct();
    model.mesh = struct();
    model.mesh.vertices = X;
    model.mesh.faces = F;
    if cindex > 1
        model.mesh.FaceVertexCData = C;
    end
end
