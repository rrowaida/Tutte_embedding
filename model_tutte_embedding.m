%
% Compute Tutte embedding of a mesh according to the uniform Laplacian
%
% function [new_model, bnd] = model_tutte_embedding(model, lap_type [, bnd])
%
% Input -
%   - model: 3D model structure
%   - lap_type: type of Laplacian to use for building the
%   parameterization. It can be one of 'uniform' or 'geometric'
%   - bnd (optional): boundary used for the parameterization, specified
%   as a list of vertex indices in the order that they appear along the
%   boundary
%
% Output -
%   - new_model: new 3D model containing and entry 'texcoord' which stores
%   the 2D coordinates of the vertices according to the parameterization
%   - bnd: list of vertices that make up the boundary used in the
%   parameterization
%
function [model, bnd, A, u] = model_tutte_embedding(model, lap_type, bnd)

    % Find a boundary of the mesh if not provided as input
    if ~exist('bnd', 'var')
        bnd = find_one_boundary(model);
        if isempty(bnd)
            disp('Could not find a boundary in this mesh!');
            return;
        end
    end

    % Map vertices along boundary to a circle of radius 0.5 centered at
    % (0.5, 0.5), so that the entire mesh fits within the range 
    % [0:1, 0:1] required for texture coordinates
    u = map_boundary_to_circle(model, bnd);

    % Create matrix for linear system
    A = param_matrix(model, lap_type, bnd);

    % Solve linear systems for u and v coordinates
    p = zeros(size(model.mesh.vertices, 1), 3);
    p(:, 1) = A \ u(:, 1);
    p(:, 2) = A \ u(:, 2);

    % Add parameterized coordinates as texture coordinates of the mesh
    model.texcoord = p;
end


% Find one possible boundary of a mesh
%
% 'bnd' stores a list of vertices in the order that they appear along
% the boundary. bnd(end) and bnd(1) connect to each other to close the
% boundary. 'bnd' is empty if no boundary can be found
%
function bnd = find_one_boundary(model)

    % Initialize the boundary
    bnd = [];

    % Find the first two vertices that are on the boundary
    vinitial = find_one_boundary_vertex(model);
    vnext = find_next_boundary_vertex(model, vinitial, -1);
    vindex = vinitial;

    % Continue traversing the boundary making sure that we always move
    % forward and not backward
    while vnext != vinitial
        fflush(stdout);
        bnd(end+1) = vnext;
        vprev = vindex;
        vindex = vnext;
        vnext = find_next_boundary_vertex(model, vindex, vprev);
        if vnext == 0
            disp('Error: unable to find next boundary vertex!');
        end
    end
end


% Find one vertex on a boundary
%
function vindex = find_one_boundary_vertex(model)
    
    % First find a boundary face
    findex = 0;
    for i = 1:size(model.mesh.faces, 1)
        if is_boundary_face(model, i)
            findex = i;
            break;
        end
    end

    % Then, find a boundary vertex in this face
    vindex = 0;
    if findex > 0
        for j = 1:3
            v1 = model.mesh.faces(findex, j);
            v2 = model.mesh.faces(findex, mod(j, 3) + 1);
            edge = [v1 v2];
            if is_boundary_edge(model, edge)
                vindex = v1;
                break;
            end
        end
    end
end


% Find a boundary vertex connected to 'vindex' which is not 'vprev', so
% that we avoid traversing back to 'vprev'
%
function vnext = find_next_boundary_vertex(model, vindex, vprev)

    vnext = 0;
    for i = 1:length(model.viv{vindex})
        v2 = model.viv{vindex}(i);
        if v2 != vprev
            edge = [vindex v2];
            if is_boundary_edge(model, edge)
                vnext = v2;
                break;
            end
        end
    end
end


% Returns true if findex is a face on a boundary
%
function t = is_boundary_face(model, findex)

    t = 0;
    if length(model.fif{findex}) < 3
        t = 1;
    end
end


% Returns true if edge is an edge on a boundary. The edge is simply a
% vector with the indices of the two vertices that form the edge
%
function t = is_boundary_edge(model, edge)

    t = 0;
    faces = intersect(model.vif{edge(1)}, model.vif{edge(2)});
    if length(faces) < 2
        t = 1;
    end
end


% Map the vertices listed in 'bnd' to a circle in the given order
%
% Map vertices along boundary to a circle of radius 0.5 centered at
% (0.5, 0.5), so that the entire mesh fits within the range [0:1, 0:1]
% required for texture coordinates
%
function u = map_boundary_to_circle(model, bnd);

    % Get length of the boundary
    l = length(bnd);

    % Initialize output coordinates
    u = zeros(size(model.mesh.vertices, 1), 2);

    % Specify points along a circle
    increment = (2*pi)/(l-1);
    for i = 1:l
        u(bnd(i), :) = [0.5*cos((i-1)*increment) + 0.5, 0.5*sin((i-1)*increment) + 0.5];
    end
end


% Create the matrix needed for Tutte embedding based on the given
% boundary
%
function A = param_matrix(model, lap_type, bnd)

    % Select the type of operator to use
    if strcmp(lap_type, 'uniform')
        A = model_uniform_laplacian(model);
    elseif strcmp(lap_type, 'geometric')
        A = model_geometric_laplacian(model);
    elseif strcmp(lap_type, 'mean')
        A = model_mean_laplacian(model);
    else
        disp(['Invalid type of Laplacian specified: ' lap_type]);
    end

    % Transform the operator into our parameterization matrix
    for i = 1:length(bnd)
        A(bnd(i), :) = 0;
        A(bnd(i), bnd(i)) = 1;
    end
end
