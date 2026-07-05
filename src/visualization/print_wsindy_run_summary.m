function print_wsindy_run_summary(details, data, cfg, mean_target_err, max_target_err, mean_full_err, max_full_err)
%PRINT_WSINDY_RUN_SUMMARY Compact one-run diagnostics for WSINDy replacement.

target_name = sprintf('Target L%d', cfg.wsindy.target_level);
test_frames = cfg.frames.test_start_idx:cfg.frames.test_end_idx;

target_norm = vecnorm(details.X_l3_true);
target_corr = details.l3_spatial_corr;
valid = isfinite(target_corr);
target_norm_valid = target_norm(valid);
target_corr_valid = target_corr(valid);

weights = target_norm_valid.^2;
if sum(weights) > 0
    weighted_target_corr = sum(weights .* target_corr_valid) / sum(weights);
else
    weighted_target_corr = NaN;
end

active = target_norm_valid > 0.05 * max(target_norm_valid);
if any(active)
    active_target_corr = mean(target_corr_valid(active));
    active_target_corr_min = min(target_corr_valid(active));
else
    active_target_corr = NaN;
    active_target_corr_min = NaN;
end

background = details.X_inside_true - details.X_l3_true;
raw_norm = vecnorm(details.X_raw_true);
mrdmd_norm = vecnorm(details.X_inside_true);
target_norm_full = vecnorm(details.X_l3_true);
background_norm = vecnorm(background);
pred_full_norm = vecnorm(background + details.X_l3_wsindy);

target_to_raw = mean(target_norm_full) / (mean(raw_norm) + eps);
background_to_raw = mean(background_norm) / (mean(raw_norm) + eps);
pred_full_to_raw = mean(pred_full_norm) / (mean(raw_norm) + eps);

mrdmd_frame_error = zeros(1, size(details.X_raw_true, 2));
for k = 1:numel(mrdmd_frame_error)
    mrdmd_frame_error(k) = norm(details.X_raw_true(:,k) - details.X_inside_true(:,k)) / ...
        (norm(details.X_raw_true(:,k)) + eps);
end

fprintf('\n================ WSINDY RUN SUMMARY ================\n');
fprintf('Dataset: %s | target: L%d | inputs: [%s]\n', ...
    cfg.data.name, cfg.wsindy.target_level, num2str(cfg.wsindy.input_levels));
fprintf('Fit frames: %d-%d | test frames: %d-%d | dt %.6g\n', ...
    cfg.frames.fit_start_idx, cfg.frames.fit_end_idx, ...
    cfg.frames.test_start_idx, cfg.frames.test_end_idx, data.dt);
fprintf('mrDMD rank: %d', cfg.mrdmd.svd_rank);
if isfield(cfg.mrdmd, 'freq_threshold_cycles_per_snapshot')
    fprintf(' | threshold %.6g cycles/snapshot\n', cfg.mrdmd.freq_threshold_cycles_per_snapshot);
elseif isfield(cfg.mrdmd, 'freq_threshold_hz')
    fprintf(' | threshold %.6g Hz\n', cfg.mrdmd.freq_threshold_hz);
else
    fprintf('\n');
end
fprintf('Selected second-pass model: %s\n', details.selected_second_pass_model);

fprintf('\nFull orchestra metrics:\n');
fprintf('  WSINDy vs raw error: mean %.2f%% | max %.2f%% (per-frame raw norm)\n', ...
    mean_full_err, max_full_err);
fprintf('  WSINDy vs raw cosine: mean %.4f | min %.4f\n', ...
    mean(details.raw_full_cos), min(details.raw_full_cos));
fprintf('  mrDMD vs raw baseline error: mean %.2f%% | max %.2f%% (per-frame raw norm)\n', ...
    mean(mrdmd_frame_error) * 100, max(mrdmd_frame_error) * 100);

fprintf('\n%s metrics:\n', target_name);
fprintf('  WSINDy vs mrDMD error: mean %.2f%% | max %.2f%% (global-scale)\n', ...
    mean_target_err, max_target_err);
fprintf('  energy-weighted spatial corr: %.4f\n', weighted_target_corr);
fprintf('  active-frame spatial corr: mean %.4f | min %.4f | active %d/%d\n', ...
    active_target_corr, active_target_corr_min, nnz(active), numel(active));
fprintf('  norm: mean %.4g | max %.4g\n', mean(target_norm_full), max(target_norm_full));

fprintf('\nCancellation balance:\n');
fprintf('  raw norm mean %.4g | mrDMD full norm mean %.4g\n', ...
    mean(raw_norm), mean(mrdmd_norm));
fprintf('  target norm mean %.4g | background norm mean %.4g | WSINDy target norm mean %.4g\n', ...
    mean(target_norm_full), mean(background_norm), mean(vecnorm(details.X_l3_wsindy)));
fprintf('  background + true target norm mean %.4g | background + WSINDy target norm mean %.4g\n', ...
    mean(vecnorm(background + details.X_l3_true)), mean(pred_full_norm));
fprintf('  target/raw %.3g | background/raw %.3g | predicted full/raw %.3g\n', ...
    target_to_raw, background_to_raw, pred_full_to_raw);

if isfield(details, 'spatial_cvx_info') && ~isempty(details.spatial_cvx_info) && ...
        isfield(details.spatial_cvx_info, 'amplitude_loss')
    fprintf('\nCVX losses:\n');
    fprintf('  derivative %.4g | amplitude %.4g | full-frame spatial %.4g\n', ...
        details.spatial_cvx_info.derivative_loss, ...
        details.spatial_cvx_info.amplitude_loss, ...
        details.spatial_cvx_info.spatial_loss);
end

fprintf('\nWorst full frames by per-frame error:\n');
if isfield(details, 'raw_full_frame_error')
    [~, order] = sort(details.raw_full_frame_error, 'descend');
    worst = order(1:min(3, numel(order)));
    for ii = 1:numel(worst)
        k = worst(ii);
        fprintf('  frame %d | err %.2f%% | cosine %.4f | target corr %.4f | raw norm %.4g | pred norm %.4g\n', ...
            test_frames(k), details.raw_full_frame_error(k) * 100, details.raw_full_cos(k), ...
            details.l3_spatial_corr(k), details.raw_frame_norm(k), details.pred_frame_norm(k));
    end
end

end
