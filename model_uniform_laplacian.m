%
% Create the uniform Laplacian operator for a mesh
%
% function L = model_uniform_laplacian(model)
%
% Input -
%   - model: 3D model structure
%
% Output -
%   - L: sparse Laplacian matrix of dimensions n x n, where n is the
%   number of vertices in the mesh
%
% See also model_geometric_laplacian
%
function L = model_uniform_laplacian(model)

    model = model_connectivity(model);
    % Get number of vertices in the mesh
    n = rows(model.mesh.vertices);

    % Initialize a sparse matrix
    L = sparse(n, n);

    % Fill matrix entries
    for i = 1:n
        % Diagonal
        L(i, i) = -1;
        % Get number of neighbors of this vertex
        num_neigh = length(model.viv{i});
        % Add entry with uniform weight for each neighbor
        L(i, model.viv{i}) = 1/num_neigh;
    end
end
