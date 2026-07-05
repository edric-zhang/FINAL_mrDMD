
close all;
clc;
clearvars -except cfg;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(project_root, 'config')));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'external', 'wsindy')));

if ~exist('cfg', 'var')
    cfg = default_experiment_config();
end

if ~isfield(cfg.wsindy, 'include_constant_term')
    cfg.wsindy.include_constant_term = false;
end
if ~isfield(cfg.wsindy, 'include_time_term')
    cfg.wsindy.include_time_term = false;
end
if ~isfield(cfg.wsindy, 'force_target_peer_terms')
    cfg.wsindy.force_target_peer_terms = false;
end
if ~isfield(cfg.wsindy, 'max_peer_l3_per_equation')
    cfg.wsindy.max_peer_l3_per_equation = 2;
end
if ~isfield(cfg.wsindy, 'second_pass_mode')
    cfg.wsindy.second_pass_mode = 'wsindy';
end
if ~isfield(cfg.wsindy, 'target_linear_ridge_alpha')
    cfg.wsindy.target_linear_ridge_alpha = 0;
end
if ~isfield(cfg.wsindy, 'target_linear_include_constant')
    cfg.wsindy.target_linear_include_constant = false;
end
if ~isfield(cfg.wsindy, 'target_linear_include_time')
    cfg.wsindy.target_linear_include_time = false;
end
if ~isfield(cfg.wsindy, 'target_linear_calibrate_amplitude')
    cfg.wsindy.target_linear_calibrate_amplitude = false;
end
if ~isfield(cfg.wsindy, 'spatial_cvx_beta')
    cfg.wsindy.spatial_cvx_beta = 1;
end
if ~isfield(cfg.wsindy, 'spatial_cvx_amp_beta')
    cfg.wsindy.spatial_cvx_amp_beta = 1;
end
if ~isfield(cfg.wsindy, 'spatial_cvx_lambda_l1')
    cfg.wsindy.spatial_cvx_lambda_l1 = 0;
end
if ~isfield(cfg.wsindy, 'spatial_cvx_sample_stride')
    cfg.wsindy.spatial_cvx_sample_stride = 5;
end

fprintf('\nActive dataset: %s\n', cfg.data.name);
fprintf('mrDMD rank: %d\n', cfg.mrdmd.svd_rank);
if isfield(cfg.mrdmd, 'freq_threshold_cycles_per_snapshot')
    fprintf('mrDMD threshold: %.6g cycles/snapshot\n', cfg.mrdmd.freq_threshold_cycles_per_snapshot);
elseif isfield(cfg.mrdmd, 'freq_threshold_hz')
    fprintf('mrDMD threshold: %.6g Hz\n', cfg.mrdmd.freq_threshold_hz);
end

data = load_spatiotemporal_dataset(cfg);

mrdmd = compute_mrdmd(data.X, data.dt, cfg);
plot_end_idx = data.m;
print_available_modes(mrdmd, cfg.frames.plot_start_idx, plot_end_idx);


%% Run one configured WSINDy model

[mean_l3_err, max_l3_err, mean_full_err, max_full_err, details] = evaluate_wsindy_parameter_set( ...
    mrdmd.list_w, mrdmd.list_b, mrdmd.list_modes, mrdmd.list_t_start, mrdmd.list_bin_widths, ...
    data.X, data.dt, data.m, data.n, cfg.wsindy.input_levels, cfg.wsindy.target_level, ...
    cfg.frames.fit_start_idx, cfg.frames.fit_end_idx, cfg.frames.test_start_idx, cfg.frames.test_end_idx, ...
    cfg.wsindy.top_input_modes_per_level, cfg.wsindy.top_target_modes, ...
    cfg.wsindy.lambda1, cfg.wsindy.lambda2, cfg.wsindy.gamma, ...
    cfg.wsindy.max_terms_per_equation, cfg.wsindy.max_quadratic_base_terms, ...
    mrdmd.list_anchor_idx, cfg.wsindy.include_constant_term, cfg.wsindy.include_time_term, ...
    cfg.wsindy.force_target_peer_terms, cfg.wsindy.max_peer_l3_per_equation, cfg.wsindy.second_pass_mode, ...
    cfg.wsindy.target_linear_ridge_alpha, cfg.wsindy.target_linear_include_constant, ...
    cfg.wsindy.target_linear_include_time, cfg.wsindy.target_linear_calibrate_amplitude, ...
    cfg.wsindy.spatial_cvx_beta, cfg.wsindy.spatial_cvx_lambda_l1, cfg.wsindy.spatial_cvx_sample_stride, ...
    cfg.wsindy.spatial_cvx_amp_beta);

prepare_result_dir(cfg.results.model_dir, cfg.results.clear_models_on_single_run);
if cfg.results.save_single_model
    stamp = char(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    model_file = fullfile(cfg.results.model_dir, ['wsindy_single_' cfg.data.name '_' stamp '.mat']);
    save(model_file, 'cfg', 'mrdmd', 'details', 'mean_l3_err', 'max_l3_err', 'mean_full_err', 'max_full_err');
end



print_wsindy_run_summary(details, data, cfg, mean_l3_err, max_l3_err, mean_full_err, max_full_err);

fprintf('\n================ SECOND-PASS TARGET EQUATIONS ================\n');
print_wsindy_equations(details.w_second, details.labels_second, details.mode_labels, details.target_cols);

plot_full_metrics(details, cfg);
plot_target_metrics(details, cfg);
plot_cancellation_metrics(details, cfg);



animate_raw_mrdmd_wsindy_comparison(details, data, cfg);
