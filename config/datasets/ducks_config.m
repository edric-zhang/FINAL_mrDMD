function cfg = ducks_config()
%DUCKS_CONFIG Scalar ducks dataset configuration.
%
% Expected data file:
%   data/processed/ducks.mat
%
% The .mat file should contain X as states x time. It can also contain x/y
% coordinates or image_shape for plotting.

cfg = default_experiment_config();

cfg = select_dataset(cfg, ...
    fullfile(cfg.project_root, 'data', 'processed', 'ducks.mat'), ...
    'scalar', ...
    'ducks');

cfg.data.plot_field = 'scalar';
cfg.data.max_working_snapshots = 400;
cfg.data.frame_start = 1;

% Start with the same model structure as CFD. Tune these after the first run.
cfg.mrdmd.L = 3;
cfg.mrdmd.freq_threshold_hz = 1000;
cfg.mrdmd.svd_rank = 25;
cfg.mrdmd.min_snapshots_per_bin = 10;

% Adjust these once you know the ducks snapshot count.
cfg.frames.fit_start_idx = 1;
cfg.frames.fit_end_idx = 100;
cfg.frames.test_start_idx = 101;
cfg.frames.test_end_idx = 150;
cfg.frames.plot_start_idx = 1;

cfg.wsindy.input_levels = [1 2];
cfg.wsindy.target_level = 3;
cfg.wsindy.top_input_modes_per_level = 5;
cfg.wsindy.top_target_modes = 10;
cfg.wsindy.lambda1 = 0.006;
cfg.wsindy.lambda2 = 0.004;
cfg.wsindy.gamma = 0.012;
cfg.wsindy.max_terms_per_equation = 6;
cfg.wsindy.max_quadratic_base_terms = 4;

end
