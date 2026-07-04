function cfg = ducks_config()
%DUCKS_CONFIG Scalar ducks dataset configuration.

cfg = default_experiment_config();

cfg = select_dataset(cfg, ...
    fullfile(cfg.project_root, 'data', 'processed', 'ducks_snapshot_matrix.mat'), ...
    'scalar', ...
    'ducks');

cfg.data.plot_field = 'scalar';
cfg.data.colormap = 'gray';
cfg.data.max_working_snapshots = 400;
cfg.data.frame_start = 1;

% From original duck video script
cfg.data.dt = 5e-4;
cfg.data.image_shape = [180 320];   % rows x cols, because reshape uses [320,180]'

cfg.mrdmd.L = 3;
cfg.mrdmd.freq_threshold_cycles_per_snapshot = 0.5;
cfg.mrdmd.svd_rank = 25;
cfg.mrdmd.min_snapshots_per_bin = 10;

cfg.frames.fit_start_idx = 301;
cfg.frames.fit_end_idx = 370;
cfg.frames.test_start_idx = 371;
cfg.frames.test_end_idx = 400;
cfg.frames.plot_start_idx = 301;

cfg.wsindy.input_levels = [1 2];
cfg.wsindy.target_level = 3;
cfg.wsindy.top_input_modes_per_level = 5;
cfg.wsindy.top_target_modes = 6;

cfg.wsindy.lambda1 = 0.002;
cfg.wsindy.lambda2 = 0.001;
cfg.wsindy.gamma = 1e-2;

cfg.wsindy.max_terms_per_equation = 8;
cfg.wsindy.max_quadratic_base_terms = 5;

end
