%
% Normalize the dimensions of a 3D model
%
% function new_model = model_normalize(model [, rng]) 
%
% Input -
%   - model: input 3D model
%   - rng (optional): target range for normalizing the shape
%
% Output -
%   - new_model is the normalized mesh. The vertex positions of the
%   longest axis are mapped to [-1.0, 1.0], while the other axes are
%   mapped so as to preserve the aspect ratio of the model.
%   Alternatively, a different target range can be provided in the 'rng'
%   parameter (a 2D vector with min and max value for the longest axis).
%
% See also model, model_read
%
function new_model = model_normalize(model, rng)
%
% Copyright (c) 2008-2018 Oliver van Kaick <ovankaic@gmail.com>
%

    % Get minimum and maximum of vertex coordinates
    mn = min(model.mesh.vertices);
    mx = max(model.mesh.vertices);

    % Define mapping parameters
    range = max(mx - mn);
    half_range = range/2;
    cnst = (mx - mn)/2;

    % Perform mapping
    if ~exist('rng', 'var')
        % Map vertex coordinates to -1.0, 1.0
        new_model = model;
        for i = 1:size(model.mesh.vertices, 1)
            new_model.mesh.vertices(i, :) = (new_model.mesh.vertices(i, :) - mn - cnst)/half_range;
        end
    else
        % Map vertex coordinates (mapping to rng)
        new_model = model;
        mn = min(mn);
        mx = max(mx);
        if length(rng(:)) == 2
            for i = 1:size(model.mesh.vertices, 1)
                new_model.mesh.vertices(i, :) = rng(1) + (rng(2) - rng(1))*((new_model.mesh.vertices(i, :) - mn)/(mx - mn));
            end
        else
            for i = 1:size(model.mesh.vertices, 1)
                for j = 1:3
                    new_model.mesh.vertices(i, j) = rng(j, 1) + (rng(j, 2) - rng(j, 1))*((new_model.mesh.vertices(i, j) - mn)/(mx - mn));
                end
            end
        end
    end
end
