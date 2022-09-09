%
% Create the geometric Laplacian operator for a mesh
%
% function L = model_geometric_laplacian(model)
%
% Input -
%   - model: 3D model structure
%
% Output -
%   - L: sparse Laplacian matrix of dimensions n x n, where n is the
%   number of vertices in the mesh
%
% See also model_uniform_laplacian
%
function [L, boundary, negative] = model_mean_laplacian(model)

    % Get number of vertices in the mesh
    model = model_connectivity(model);
    n = rows(model.mesh.vertices);
    
    % Initialize output
    L = sparse(n, n);

    % Fill matrix entries
    for i = 1:n % for all vertices
        
        x_i_j = 0;
        %sum = 0;
        % Get vertex coordinates
        v = model.mesh.vertices(i, :);

        % For each neighboring vertex
        for j = 1:length(model.viv{i})
          
            % total number of neighbors
            
            n_num=length(model.viv{i});
            
            % Get neighboring vertex
            % (i, neigh) is the edge we are considering
            neigh_index = model.viv{i}(j);
            neigh = model.mesh.vertices(neigh_index, :);
            
            % Get two other vertices that are not i or neigh
            faces = intersect(model.vif{i}, model.vif{neigh_index});
            if length(faces) ~= 2
                boundary = 1;
            end
            ov = zeros(0, 3);
            for k = 1:length(faces)
                for l = 1:3
                    if (model.mesh.faces(faces(k), l) != i) && ...
                       (model.mesh.faces(faces(k), l) != neigh_index)
                       ov(end+1, :) = model.mesh.vertices(model.mesh.faces(faces(k), l), :);
                    end
                end
            end

            % Now compute cotangents for mean curvature and the estimation
            % of Voronoi area
            cot = [0 0];
            for k = 1:size(ov, 1)
                % Compute the cotangent of the corner k
                vec1 = ov(k, :) - v;
                vec1 = vec1 ./ norm(vec1, 2);
                vec2 = neigh - v;
                vec2 = vec2 ./ norm(vec2, 2);
                cosine = dot(vec1, vec2);
                sine = norm(cross(vec1, vec2), 2);
                s = sin(asin(sine)/2);
                c = cos(acos(cosine)/2);
                cot(k) = s/c;
                % Check if we have an obtuse angle
                if cot(k) < 0
                    % Use mixed area element
                    % Compute area of triangle
                    % cross / 2 to get area of triangle from area of parallelogram
                    % result / 2 to get 1/2 of the triangle area
                    cot(k) = norm(cross(ov(k, :) - v, neigh - v), 2)/4;
                end
            end
            L(i, model.viv{i}(j)) = cot(1) + cot(2);
            x_i_j = sqrt((v - neigh)*(v - neigh)');
            
            % Check negativity of cotan weights
            if (cot(1) + cot(2)) < 0
                negative = 1;
            end
        end

        % Weight all entries by area
        L(i, :) = (1/x_i_j) * L(i, :);

        % Diagonal
        L(i, i) = -sum(L(i, :));
    end
end
