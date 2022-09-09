% Load and prepare the mesh
model = model_read('meshes/bunny.obj');
model = model_normalize(model);
model = model_connectivity(model);

% Perform smoothing
lap_type = {'uniform', 'geometric'};
for i = 1:length(lap_type)
    % Build operator
    if i == 1
        L = model_uniform_laplacian(model);
    else
        L = model_geometric_laplacian(model);
    end

    % Apply operator
    num_iterations = 20;
    for n = 1:num_iterations
        model.mesh.vertices = model.mesh.vertices + L*model.mesh.vertices;
    end

    % Save mesh
    model_write(model, ['bunny_smooth_' lap_type{i} '.obj']);
end
