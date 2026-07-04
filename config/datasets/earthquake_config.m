function cfg = earthquake_config()
%EARTHQUAKE_CONFIG Scalar earthquake dataset configuration.
%
% Expected data file:
%   data/processed/earthquake.mat
%
% The .mat file should contain X as states x time. It can also contain x/y
% coordinates or image_shape for plotting.

cfg = default_experiment_config();

cfg = select_dataset(cfg, ...
    fullfile(cfg.project_root, 'data', 'processed', 'earthquake.mat'), ...
    'scalar', ...
    'earthquake');

cfg.data.plot_field = 'scalar';
cfg.data.colormap = 'hot';
cfg.data.nan_fill_value = 0;
cfg.data.max_working_snapshots = 251;
cfg.data.frame_start = 1;

cfg.mrdmd.L = 5;
cfg.mrdmd.freq_threshold_cycles_per_snapshot = 0.16;

cfg.mrdmd.svd_rank = 75;
cfg.mrdmd.min_snapshots_per_bin = 10;

% mrDMD can use more than 3 levels, but WSINDy is still trained on level 3
% because those bins have enough snapshots for fit/test.
% For 251 snapshots, the level-3 bins are approximately:
%   L3 B1: 1-62, B2: 63-125, B3: 126-188, B4: 189-251.
% Keep fit/test inside one level-3 bin for the WSINDy recreation.
cfg.frames.fit_start_idx = 189;
cfg.frames.fit_end_idx = 230;
cfg.frames.test_start_idx = 231;
cfg.frames.test_end_idx = 251;
cfg.frames.plot_start_idx = 189;

cfg.wsindy.input_levels = [1 2];
cfg.wsindy.target_level = 3;
cfg.wsindy.top_input_modes_per_level = 5;
cfg.wsindy.top_target_modes = 1;
cfg.wsindy.lambda1 = 0.4;
cfg.wsindy.lambda2 = 0.08;
cfg.wsindy.gamma = 0.001;
cfg.wsindy.max_terms_per_equation = 5;
cfg.wsindy.max_quadratic_base_terms = 4;

cfg.results.save_sweep_results = false;

end
