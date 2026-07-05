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
data = load_spatiotemporal_dataset(cfg);

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

fprintf('\nSTARTING MRDMD\n');
mrdmd = compute_mrdmd(data.X, data.dt, cfg);

[~, mrdmd_error] = reconstruct_mrdmd(mrdmd, data.X, cfg.data.frame_start, data.numsnapshots);
plot_reconstruction_error(mrdmd_error, data.numsnapshots);

plot_end_idx = data.m;
print_available_modes(mrdmd, cfg.frames.plot_start_idx, plot_end_idx);

% Format: [Level, Bin (j), Mode_Index] (Use 0 as a wildcard for ALL)
% Use 1,0,0 - 2,0,0 - 3,0,0 for full period-specific reconstruction
target_coordinates = [
    %1, 0, 0;
    %2, 0, 0;
    %3, 0, 0;
    3, 0, 0;
];

%% Run WSINDy structure/parameter sweep

fprintf('\nRunning weak SINDy structure sweep on mrDMD modal amplitudes...\n');
tic;

param_grid = sweep_grid_cfd();
sweep_results = zeros(size(param_grid,1), 15);

for rr = 1:size(param_grid,1)
    lambda1 = param_grid(rr,1);
    lambda2 = param_grid(rr,2);
    gamma = param_grid(rr,3);
    top_input_modes_per_level = param_grid(rr,4);
    top_target_modes = param_grid(rr,5);
    max_terms_per_equation = param_grid(rr,6);
    max_quadratic_base_terms = param_grid(rr,7);

    fprintf('\n=== RUN %d/%d | lambda1 %.4g | lambda2 %.4g | gamma %.4g | inputs %d | targets %d | maxTerms %d | quadBase %d ===\n', ...
        rr, size(param_grid,1), lambda1, lambda2, gamma, ...
        top_input_modes_per_level, top_target_modes, max_terms_per_equation, max_quadratic_base_terms);

    [mean_l3_err, max_l3_err, mean_full_err, max_full_err, details_sweep] = evaluate_wsindy_parameter_set( ...
        mrdmd.list_w, mrdmd.list_b, mrdmd.list_modes, mrdmd.list_t_start, mrdmd.list_bin_widths, ...
        data.X, data.dt, data.m, data.n, cfg.wsindy.input_levels, cfg.wsindy.target_level, ...
        cfg.frames.fit_start_idx, cfg.frames.fit_end_idx, cfg.frames.test_start_idx, cfg.frames.test_end_idx, ...
        top_input_modes_per_level, top_target_modes, ...
        lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms, ...
        mrdmd.list_anchor_idx, cfg.wsindy.include_constant_term, cfg.wsindy.include_time_term, ...
    cfg.wsindy.force_target_peer_terms, cfg.wsindy.max_peer_l3_per_equation, cfg.wsindy.second_pass_mode, ...
    cfg.wsindy.target_linear_ridge_alpha, cfg.wsindy.target_linear_include_constant, ...
    cfg.wsindy.target_linear_include_time, cfg.wsindy.target_linear_calibrate_amplitude, ...
    cfg.wsindy.spatial_cvx_beta, cfg.wsindy.spatial_cvx_lambda_l1, cfg.wsindy.spatial_cvx_sample_stride, ...
    cfg.wsindy.spatial_cvx_amp_beta);

    mean_full_corr = mean(details_sweep.raw_full_corr);
    min_full_corr = min(details_sweep.raw_full_corr);
    mean_modal_corr = mean(details_sweep.mrdmd_full_corr);
    mean_mrdmd_raw_err = mean(details_sweep.mrdmd_raw_error) * 100;

    sweep_results(rr,:) = [lambda1, lambda2, gamma, ...
        top_input_modes_per_level, top_target_modes, max_terms_per_equation, max_quadratic_base_terms, ...
        mean_l3_err, max_l3_err, mean_full_err, max_full_err, ...
        mean_full_corr, min_full_corr, mean_modal_corr, mean_mrdmd_raw_err];

    fprintf('RESULT | lambda1 %.4g | lambda2 %.4g | gamma %.4g | inputs %d | targets %d | maxTerms %d | quadBase %d | L3 mean %.2f%% | L3 max %.2f%% | raw full mean %.2f%% | raw full max %.2f%% | raw corr mean %.4f | raw corr min %.4f\n', ...
        lambda1, lambda2, gamma, top_input_modes_per_level, top_target_modes, ...
        max_terms_per_equation, max_quadratic_base_terms, ...
        mean_l3_err, max_l3_err, mean_full_err, max_full_err, mean_full_corr, min_full_corr);
end

fprintf('\n================ STRUCTURE SWEEP RESULTS ================\n');
result_table = array2table(sweep_results, ...
    'VariableNames', {'lambda1','lambda2','gamma','top_inputs','top_targets','max_terms','quad_base', ...
    'mean_L3_error','max_L3_error','mean_raw_full_error','max_raw_full_error', ...
    'mean_raw_full_corr','min_raw_full_corr','mean_mrdmd_full_corr','mean_mrdmd_raw_error'});
disp(result_table);

prepare_result_dir(cfg.results.sweep_dir, cfg.results.clear_sweeps_on_sweep_run);
if cfg.results.save_sweep_results
    stamp = char(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    mat_file = fullfile(cfg.results.sweep_dir, ['wsindy_sweep_' cfg.data.name '_' stamp '.mat']);
    csv_file = fullfile(cfg.results.sweep_dir, ['wsindy_sweep_' cfg.data.name '_' stamp '.csv']);
    save(mat_file, 'cfg', 'sweep_results', 'result_table');
    writetable(result_table, csv_file);

    fprintf('\nSaved sweep MAT: %s\n', mat_file);
    fprintf('Saved sweep CSV: %s\n', csv_file);
else
    fprintf('\nSweep saving is off. Set cfg.results.save_sweep_results = true to save sweep files.\n');
end
toc;
