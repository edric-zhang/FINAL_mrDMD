function cfg = default_experiment_config()
%DEFAULT_EXPERIMENT_CONFIG Baseline settings for CFD mrDMD + WSINDy runs.

cfg.project_root = fileparts(fileparts(mfilename('fullpath')));

cfg.data.mat_file = fullfile(cfg.project_root, 'data', 'processed', 'cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat');
cfg.data.name = 'cylinder_uv';
% State layout options:
%   'uv'     means X is stacked as [u; v], so size(X,1) = 2*npoints.
%   'scalar' means X has one scalar state per spatial point/pixel.
%   'auto'   tries 'uv' when npoints is supplied and size(X,1)=2*npoints,
%            otherwise treats X as scalar.
cfg.data.state_layout = 'uv';
cfg.data.plot_field = 'speed'; % For 'uv': 'u', 'v', or 'speed'. For 'scalar': 'scalar'.
cfg.data.max_working_snapshots = 400;
cfg.data.frame_start = 1;

cfg.mrdmd.L = 3;
cfg.mrdmd.freq_threshold_cycles_per_snapshot = 0.5;
cfg.mrdmd.svd_rank = 25;
cfg.mrdmd.min_snapshots_per_bin = 10;

% MULTI-TARGET EXTRACTION WITH ACTIVE LOOKUP MENU
cfg.frames.fit_start_idx = 263;
cfg.frames.fit_end_idx = 310;
cfg.frames.test_start_idx = 311;
cfg.frames.test_end_idx = 350;
cfg.frames.plot_start_idx = 100;

% WEAK SINDY ON MRDMD TEMPORAL COEFFICIENTS
% Goal:
% Learn d/dt of Level 3 modal amplitudes using Level 1-3 modal amplitudes.
%
% This replaces the old spatial regression:
%   L3 spatial modes = Theta(L1-L2 spatial modes)*Xi
%
% with actual dynamics:
%   d(a_L3)/dt = Theta(a_L1,a_L2,a_L3)*Xi
cfg.wsindy.input_levels = [1 2];
cfg.wsindy.target_level = 3;
cfg.wsindy.top_input_modes_per_level = 5;
cfg.wsindy.top_target_modes = 10;
cfg.wsindy.lambda1 = 0.006;
cfg.wsindy.lambda2 = 0.004;
cfg.wsindy.gamma = 0.012;
cfg.wsindy.max_terms_per_equation = 6;
cfg.wsindy.max_quadratic_base_terms = 4;
cfg.wsindy.include_constant_term = false;
cfg.wsindy.include_time_term = false;
cfg.wsindy.force_target_peer_terms = false;
cfg.wsindy.max_peer_l3_per_equation = 2;
cfg.wsindy.second_pass_mode = 'wsindy';
cfg.wsindy.target_linear_ridge_alpha = 0;
cfg.wsindy.target_linear_include_constant = false;
cfg.wsindy.target_linear_include_time = false;
cfg.wsindy.target_linear_calibrate_amplitude = false;
cfg.wsindy.spatial_cvx_beta = 1;
cfg.wsindy.spatial_cvx_amp_beta = 1;
cfg.wsindy.spatial_cvx_lambda_l1 = 0;
cfg.wsindy.spatial_cvx_sample_stride = 5;
cfg.results.sweep_dir = fullfile(cfg.project_root, 'results', 'sweeps');
cfg.results.figure_dir = fullfile(cfg.project_root, 'results', 'figures');
cfg.results.model_dir = fullfile(cfg.project_root, 'results', 'models');
cfg.results.save_single_model = true;
cfg.results.save_sweep_results = false;
% Turn these on only when you intentionally want to clear old generated files.
cfg.results.clear_models_on_single_run = false;
cfg.results.clear_sweeps_on_sweep_run = false;

end
