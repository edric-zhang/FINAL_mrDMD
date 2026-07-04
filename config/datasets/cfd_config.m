function cfg = cfd_config()
%CFD_CONFIG Current cylinder wake [u; v] CFD dataset configuration.

cfg = default_experiment_config();

cfg = select_dataset(cfg, ...
    fullfile(cfg.project_root, 'data', 'processed', 'cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat'), ...
    'uv', ...
    'cylinder_uv');

% For this dataset, use 'speed', 'u', or 'v'.
cfg.data.plot_field = 'speed';
cfg.data.max_working_snapshots = 400;
cfg.data.frame_start = 1;

cfg.mrdmd.L = 3;
cfg.mrdmd.freq_threshold_cycles_per_snapshot = 0.5;
cfg.mrdmd.svd_rank = 25;
cfg.mrdmd.min_snapshots_per_bin = 10;

cfg.frames.fit_start_idx = 263;
cfg.frames.fit_end_idx = 310;
cfg.frames.test_start_idx = 311;
cfg.frames.test_end_idx = 350;
cfg.frames.plot_start_idx = 100;

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
