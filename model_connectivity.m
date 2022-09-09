%
% Obtain mesh connectivity information from a 3D model
%
% function new_model = model_connectivity(model);
%
% Input -
%   - model: input 3D model
%
% Output -
%   - new_model: new 3D model containing entries 'vif', 'fif', and
%   'viv'. These entries are as follows:
%       = vif: faces incident to vertices. The entry vif{i} is a vector
%       containing the indices of faces that are incident to vertex 'i'
%       (that is, faces that contain vertex 'i') -- this is the relation V -> F
%       = fif: faces adjacent to faces. The entry fif{i} is a vector
%       containing the indices of faces that are neighbors of face 'i' (that
%       is, faces which have two vertices in common with face 'i') -- this
%       is the relation F -> F
%       = viv: vertices adjacent to vertices. The entry viv{i} is a vector
%       containing the indices of vertices that are adjacent to vertex 'i'
%       (that is, vertices that have a face in common with vertex 'i') --
%       this is the relation V -> V
%       - The relation F -> V is given by the mesh itself
%
% See also model
%
function new_model = model_connectivity(model);
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>
%

    % Initialize output
    new_model = model;

    % Call the respective functions for each neighborhood relation
    new_model.vif = get_v_inc_f(model.mesh.vertices, model.mesh.faces);
    new_model.fif = get_f_adj_f(model.mesh.vertices, model.mesh.faces, new_model.vif);
    new_model.viv = get_v_adj_v(model.mesh.vertices, model.mesh.faces, new_model.vif);
end


%
% This function obtains all the faces that are incident to the vertices
% of the mesh.
%
% Input -
%   - V: mesh vertices: V(i, :) represents the 3D coordinates of vertex
%   'i'
%   - F: mesh faces: F(i, :) contains the indices of the three vertices
%   that compose face 'i'
%
% Output -
%   - VIF: incident faces: VIF{i} is a vector containing the indices of
%   faces that are incident to vertex 'i' (that is, faces that contain
%   vertex 'i')
%
function VIF = get_v_inc_f(V, F)

    % Allocate cell complex
    VIF = cell(size(V, 1), 1);

    % For each face
    for i = 1:size(F, 1)
        % For each of its vertices
        for j = 1:3
            % Add face 'i' to the list of faces of vertex F(i, j)
            VIF{F(i, j)}(end+1, 1) = i;
        end
    end
end


%
% This function obtains all the adjacencies between faces of the mesh
% (neighboring faces).
%
% Input -
%   - V: mesh vertices: V(i, :) represents the 3D coordinates of vertex
%   'i'
%   - F: mesh faces: F(i, :) contains the indices of the three vertices
%   that compose face 'i'
%   - VIF: adjacent faces: VIF{i} is a vector containing the indices of
%   faces that are incident to vertex 'i', as returned by function
%   get_v_inc_f()
%
% Output -
%   - FIF: adjacent faces: FIF{i} is a vector containing the indices of
%   faces that are neighbors of face 'i' (that is, faces which have two
%   vertices in common with face 'i')
%
function FIF = get_f_adj_f(V, F, VIF)

    % Init cell complex
    FIF = cell(size(F, 1), 1);

    % For each face
    for i = 1:size(F, 1)
        % For each vertex of the face
        for j = 1:3
            % For each face incident to the vertex
            for k = 1:size(VIF{F(i, j)}, 1)
                % If the face is not itself
                if VIF{F(i, j)}(k) ~= i
                    % If it is in the list of another vertex
                    other1 = mod(j, 3) + 1;
                    other2 = mod(j + 1, 3) + 1;
                    if (sum(VIF{F(i, j)}(k) == VIF{F(i, other1)}) > 0) || ...
                       (sum(VIF{F(i, j)}(k) == VIF{F(i, other2)}) > 0) 
                        % If it was not already added to the vector
                        if sum(FIF{i} == VIF{F(i, j)}(k)) == 0
                            % Add face to the end of the vector
                            FIF{i}(end+1, 1) = VIF{F(i, j)}(k);
                        end
                    end
                end
            end
        end
    end
end


%
% This function obtains all the adjacencies between vertices of the mesh
%
% Input -
%   - V: mesh vertices: V(i, :) represents the 3D coordinates of vertex
%   'i'
%   - F: mesh faces: F(i, :) contains the indices of the three vertices
%   that compose face 'i'
%   - VIF: incident faces: VIF{i} is a vector containing the indices of
%   faces that are incident to vertex 'i', as returned by function
%   get_v_inc_f()
%
% Output -
%   - VIV: adjacent vertices: VIV{i} is a vector containing the indices
%   of vertices that are adjacent to vertex 'i' (that is, vertices that
%   have a face in common with vertex 'i')
%
function VIV = get_v_adj_v(V, F, VIF)

    % Init cell complex
    VIV = cell(size(V, 1), 1);

    % For each vertex
    for i = 1:size(V, 1)
        % For each of its incident faces
        for j = 1:size(VIF{i}, 1)
            % For each vertex of the incident face
            for k = 1:3
                % If vertex is not itself
                if F(VIF{i}(j), k) ~= i
                    % If it was not already added to the vector
                    if sum(VIV{i} == F(VIF{i}(j), k)) == 0
                        % Add vertex to the end of the vector
                        VIV{i}(end+1, 1) = F(VIF{i}(j), k);
                    end
                end
            end
        end
    end
end
