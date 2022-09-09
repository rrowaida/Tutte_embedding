% Specify input and output
%input_filename = 'meshes/crater.obj';
%output_prefix_texture = 'crater_textured';
%output_prefix_param = 'crater_param';
%output_prefix_boundary = 'crater_boundary';
%input_filename = 'camel_cut.obj';
%output_prefix_texture = 'camel_textured';
%output_prefix_param = 'camel_param';
%output_prefix_boundary = 'camel_boundary';
input_filename = 'beetle.obj';
output_prefix_texture = 'beetle_textured';
output_prefix_param = 'beetle_param';
output_prefix_boundary = 'beetle_boundary';

% Load input model if needed
if 1
model = model_read(input_filename);
model = model_normalize(model);
model = model_connectivity(model);
end

% Compute parameterizations with different operators
lap_type = {'uniform', 'geometric', 'mean'};

for i = 1:length(lap_type)

    % Compute Tutte embedding
    [model, bnd, A, u] = model_tutte_embedding(model, lap_type{i});

    % Transform the boundary into colors, if we would like to visualize the
    % boundary on the mesh
    model.mesh.FaceVertexCData = ones(size(model.mesh.vertices, 1), 3);
    model.mesh.FaceVertexCData(bnd, :) = repmat([1 0 0], length(bnd), 1);
    model_write(model, [output_prefix_boundary '_' lap_type{i} '.obj']);
    model.mesh = rmfield(model.mesh, 'FaceVertexCData');

    % Save texture mesh as output
    model_write_obj_texcoord(model, [output_prefix_texture '_' lap_type{i} '.obj']);

    % Save mesh parametrized to 2D
    param_model = model;
    param_model.mesh.vertices = param_model.texcoord;
    param_model = rmfield(param_model, 'texcoord');
    model_write(param_model, [output_prefix_param '_' lap_type{i} '.obj']);
end
