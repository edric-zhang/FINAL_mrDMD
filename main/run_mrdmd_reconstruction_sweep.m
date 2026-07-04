close all;
clc;
clearvars -except cfg mrdmd_grid;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(project_root, 'config')));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'external', 'wsindy')));

if ~exist('cfg', 'var')
    cfg = earthquake_config();
end

if ~exist('mrdmd_grid', 'var') || isempty(mrdmd_grid)
    % rank, L, freq_threshold_cycles_per_snapshot
    mrdmd_grid = [
        50   4   0.08
        50   4   0.12
        50   4   0.16
        50   5   0.08
        50   5   0.12
        50   5   0.16

        75   4   0.08
        75   4   0.12
        75   4   0.16
        75   5   0.08
        75   5   0.12
        75   5   0.16

        100  4   0.08
        100  4   0.12
        100  4   0.16
        100  5   0.08
        100  5   0.12
        100  5   0.16
    ];
end

data = load_spatiotemporal_dataset(cfg);
test_idx = cfg.frames.test_start_idx:cfg.frames.test_end_idx;

fprintf('\n================ MRDMD RECONSTRUCTION SWEEP ================\n');
fprintf('Dataset: %s | Test frames: %d-%d\n', cfg.data.name, test_idx(1), test_idx(end));

num_runs = size(mrdmd_grid, 1);
sweep_results = zeros(num_runs, 9);

for rr = 1:num_runs
    rank_val = mrdmd_grid(rr, 1);
    L_val = mrdmd_grid(rr, 2);
    freq_val = mrdmd_grid(rr, 3);

    cfg_run = cfg;
    cfg_run.mrdmd.svd_rank = rank_val;
    cfg_run.mrdmd.L = L_val;
    cfg_run.mrdmd.freq_threshold_cycles_per_snapshot = freq_val;

    fprintf('\n=== MRDMD RUN %d/%d | rank %d | L %d | freq %.4g ===\n', ...
        rr, num_runs, rank_val, L_val, freq_val);

    tic;
    mrdmd = compute_mrdmd(data.X, data.dt, cfg_run);
    [X_rec, mrdmd_error] = reconstruct_mrdmd(mrdmd, data.X, ...
        cfg_run.data.frame_start, data.numsnapshots);
    elapsed_sec = toc;

    corr_vals = zeros(1, numel(test_idx));
    rel_corr_vals = zeros(1, numel(test_idx));

    for kk = 1:numel(test_idx)
        frame = test_idx(kk);
        raw = data.X(:, frame);
        rec = X_rec(:, frame);
        raw_centered = raw - mean(raw);
        rec_centered = rec - mean(rec);

        corr_vals(kk) = dot(raw, rec) / ((norm(raw) * norm(rec)) + eps);
        rel_corr_vals(kk) = dot(raw_centered, rec_centered) / ...
            ((norm(raw_centered) * norm(rec_centered)) + eps);
    end

    test_error_pct = mrdmd_error(test_idx) * 100;
    all_error_pct = mrdmd_error(~isnan(mrdmd_error)) * 100;

    mean_test_corr = mean(corr_vals);
    min_test_corr = min(corr_vals);
    mean_centered_corr = mean(rel_corr_vals);
    mean_test_error = mean(test_error_pct);
    max_test_error = max(test_error_pct);
    mean_all_error = mean(all_error_pct);

    sweep_results(rr,:) = [rank_val, L_val, freq_val, ...
        mean_test_corr, min_test_corr, mean_centered_corr, ...
        mean_test_error, max_test_error, mean_all_error];

    fprintf(['RESULT | rank %d | L %d | freq %.4g | mean corr %.4f | min corr %.4f | ' ...
        'centered corr %.4f | test err %.2f%% | max test err %.2f%% | all err %.2f%% | %.2fs\n'], ...
        rank_val, L_val, freq_val, mean_test_corr, min_test_corr, mean_centered_corr, ...
        mean_test_error, max_test_error, mean_all_error, elapsed_sec);
end

result_table = array2table(sweep_results, ...
    'VariableNames', {'rank','L','freq_threshold','mean_test_corr','min_test_corr', ...
    'mean_centered_corr','mean_test_error','max_test_error','mean_all_error'});

result_table.score = result_table.mean_test_error - 100 * result_table.mean_test_corr;
result_table = sortrows(result_table, 'score');

fprintf('\n================ MRDMD RECONSTRUCTION SWEEP RESULTS ================\n');
disp(result_table);

prepare_result_dir(cfg.results.sweep_dir, cfg.results.clear_sweeps_on_sweep_run);
if isfield(cfg.results, 'save_sweep_results') && cfg.results.save_sweep_results
    stamp = char(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    mat_file = fullfile(cfg.results.sweep_dir, ['mrdmd_reconstruction_sweep_' cfg.data.name '_' stamp '.mat']);
    csv_file = fullfile(cfg.results.sweep_dir, ['mrdmd_reconstruction_sweep_' cfg.data.name '_' stamp '.csv']);
    save(mat_file, 'cfg', 'mrdmd_grid', 'sweep_results', 'result_table');
    writetable(result_table, csv_file);
    fprintf('\nSaved sweep MAT: %s\n', mat_file);
    fprintf('Saved sweep CSV: %s\n', csv_file);
end
