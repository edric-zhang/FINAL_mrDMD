function cfg = select_dataset(cfg, mat_file, state_layout, dataset_name)
%SELECT_DATASET Point a config at a different spatiotemporal dataset.
%
% Example:
%   cfg = default_experiment_config();
%   cfg = select_dataset(cfg, ...
%       fullfile(cfg.project_root, 'data', 'processed', 'my_scalar_data.mat'), ...
%       'scalar', 'my_scalar_data');

cfg.data.mat_file = mat_file;
cfg.data.state_layout = state_layout;

switch lower(state_layout)
    case 'uv'
        cfg.data.plot_field = 'speed';
    case 'scalar'
        cfg.data.plot_field = 'scalar';
    otherwise
        cfg.data.plot_field = 'speed';
end

if nargin >= 4 && ~isempty(dataset_name)
    cfg.data.name = dataset_name;
else
    [~, cfg.data.name] = fileparts(mat_file);
end

end
