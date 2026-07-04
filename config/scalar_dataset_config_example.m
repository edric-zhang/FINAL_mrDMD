function cfg = scalar_dataset_config_example(mat_file)
%SCALAR_DATASET_CONFIG_EXAMPLE Example config for one-state-per-pixel data.
%
% Copy this pattern when you want to run a non-[u;v] dataset.

cfg = default_experiment_config();
cfg = select_dataset(cfg, mat_file, 'scalar');

% Tune these for the scalar dataset length.
cfg.data.max_working_snapshots = 400;
cfg.frames.fit_start_idx = 263;
cfg.frames.fit_end_idx = 310;
cfg.frames.test_start_idx = 311;
cfg.frames.test_end_idx = 350;
cfg.frames.plot_start_idx = 100;

end
