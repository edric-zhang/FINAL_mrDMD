function [mean_l3_err, max_l3_err, mean_full_err, max_full_err, details] = evaluate_wsindy_parameter_set( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    X, dt, m, n, input_levels, target_level, ...
    fit_start_idx, fit_end_idx, test_start_idx, test_end_idx, ...
    top_input_modes_per_level, top_target_modes, ...
    lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms, ...
    list_anchor_idx, include_constant_term, include_time_term, force_target_peer_terms, max_peer_l3_per_equation, second_pass_mode, ...
    target_linear_ridge_alpha, target_linear_include_constant, target_linear_include_time, target_linear_calibrate_amplitude, ...
    spatial_cvx_beta, spatial_cvx_lambda_l1, spatial_cvx_sample_stride, ...
    spatial_cvx_amp_beta)
%EVALUATE_WSINDY_PARAMETER_SET Train and score one WSINDy parameter set.

if nargin < 24 || isempty(include_constant_term)
    include_constant_term = false;
end
if nargin < 25 || isempty(include_time_term)
    include_time_term = false;
end
if nargin < 26 || isempty(force_target_peer_terms)
    force_target_peer_terms = false;
end
if nargin < 27 || isempty(max_peer_l3_per_equation)
    max_peer_l3_per_equation = 2;
end
if nargin < 28 || isempty(second_pass_mode)
    second_pass_mode = 'wsindy';
end
if nargin < 29 || isempty(target_linear_ridge_alpha)
    target_linear_ridge_alpha = 0;
end
if nargin < 30 || isempty(target_linear_include_constant)
    target_linear_include_constant = false;
end
if nargin < 31 || isempty(target_linear_include_time)
    target_linear_include_time = false;
end
if nargin < 32 || isempty(target_linear_calibrate_amplitude)
    target_linear_calibrate_amplitude = false;
end
if nargin < 33 || isempty(spatial_cvx_beta)
    spatial_cvx_beta = 1;
end
if nargin < 34 || isempty(spatial_cvx_lambda_l1)
    spatial_cvx_lambda_l1 = 0;
end
if nargin < 35 || isempty(spatial_cvx_sample_stride)
    spatial_cvx_sample_stride = 5;
end
if nargin < 36 || isempty(spatial_cvx_amp_beta)
    spatial_cvx_amp_beta = 1;
end

[w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs] = run_mrdmd_wsindy( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    dt, m, input_levels, target_level, fit_start_idx, fit_end_idx, ...
    top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
    max_terms_per_equation, max_quadratic_base_terms, list_anchor_idx, include_constant_term, include_time_term, force_target_peer_terms, max_peer_l3_per_equation);

test_steps = test_end_idx - test_start_idx + 1;
recreate_steps = test_steps;

[xobs_test, tobs_test] = build_mrdmd_xobs_for_labels( ...
    list_w, list_b, list_t_start, list_bin_widths, ...
    dt, m, mode_labels, test_start_idx, test_end_idx, list_anchor_idx);

xobs_test = (xobs_test - mu_xobs) ./ sigma_xobs;

% One-step teacher-forced derivative validation.
% Use the true current modal state at frame q and ask the learned
% derivative to predict only frame q+1. This directly tests whether the
% derivative equations produce visible local motion before free-running.
y_recreate = xobs_test(:, target_cols);

for q = 1:(test_steps-1)
    y_now = xobs_test(q, target_cols).';

    dydt_now = mrdmd_wsindy_inside_rhs( ...
        tobs_test(q), y_now, ...
        w_second, labels_second, mode_labels, target_cols, ...
        xobs_test, tobs_test);

    y_recreate(q+1, :) = xobs_test(q, target_cols) + dt * dydt_now.';
end

y_recreate_phys = y_recreate .* sigma_xobs(target_cols) + mu_xobs(target_cols);

X_l3_true = build_true_l3_contribution_from_modes( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    list_anchor_idx, mode_labels, target_cols, test_start_idx, recreate_steps, n, X);

snapshot_norms = vecnorm(X);
valid_norms = snapshot_norms(snapshot_norms > 0);
global_scale = median(valid_norms);

% Full mrDMD background is needed before candidate selection because a good
% isolated L3 match can still destroy the raw frame when L3 cancels other levels.
mrdmd_tmp.list_w = list_w;
mrdmd_tmp.list_b = list_b;
mrdmd_tmp.list_modes = list_modes;
mrdmd_tmp.list_t_start = list_t_start;
mrdmd_tmp.list_bin_widths = list_bin_widths;
mrdmd_tmp.list_anchor_idx = list_anchor_idx;
mrdmd_tmp.L = size(list_modes, 1);

[X_mrdmd_full, ~] = reconstruct_mrdmd(mrdmd_tmp, X, 1, size(X, 2));
fit_steps = fit_end_idx - fit_start_idx + 1;
X_l3_fit_true = build_true_l3_contribution_from_modes( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    list_anchor_idx, mode_labels, target_cols, fit_start_idx, fit_steps, n, X);
X_fit_background = real(X_mrdmd_full(:, fit_start_idx:fit_end_idx)) - real(X_l3_fit_true);
X_raw_fit = real(X(:, fit_start_idx:fit_end_idx));

candidate_report = struct('name', {}, 'mean_l3_corr', {}, 'mean_l3_error', {}, 'max_l3_error', {});

[X_l3_wsindy, l3_recreate_error, l3_corr] = score_l3_candidate( ...
    y_recreate_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X);

best_name = 'wsindy';
best_y_recreate_phys = y_recreate_phys;
best_X_l3_wsindy = X_l3_wsindy;
best_l3_recreate_error = l3_recreate_error;
best_l3_corr = l3_corr;
best_w_second = w_second;
best_labels_second = labels_second;
best_calibration_report = struct('label', {}, 'applied', {}, 'fit_amp_corr', {}, 'fit_rel_err', {}, 'scale', {}, 'bias', {});

candidate_report(end+1) = make_candidate_report(best_name, best_l3_corr, best_l3_recreate_error);
spatial_cvx_info = [];

if strcmpi(second_pass_mode, 'spatial_cvx')
    [xobs_fit_raw, tobs_fit] = build_mrdmd_xobs_for_labels( ...
        list_w, list_b, list_t_start, list_bin_widths, ...
        dt, m, mode_labels, fit_start_idx, fit_end_idx, list_anchor_idx);

    xobs_fit = (xobs_fit_raw - mu_xobs) ./ sigma_xobs;

    [w_spatial, spatial_cvx_info] = run_mrdmd_wsindy_spatial_cvx( ...
        labels_second, mode_labels, target_cols, ...
        xobs_fit, tobs_fit, mu_xobs, sigma_xobs, list_modes, ...
        spatial_cvx_beta, spatial_cvx_lambda_l1, spatial_cvx_sample_stride);

    y_spatial = predict_wsindy_one_step(xobs_test, tobs_test, target_cols, ...
        w_spatial, labels_second, mode_labels, dt);
    y_spatial_phys = y_spatial .* sigma_xobs(target_cols) + mu_xobs(target_cols);

    [X_spatial, err_spatial, corr_spatial] = score_l3_candidate( ...
        y_spatial_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X);

    candidate_report(end+1) = make_candidate_report('spatial_cvx', corr_spatial, err_spatial);

    best_name = 'spatial_cvx';
    best_y_recreate_phys = y_spatial_phys;
    best_X_l3_wsindy = X_spatial;
    best_l3_recreate_error = err_spatial;
    best_l3_corr = corr_spatial;
    best_w_second = w_spatial;
    best_labels_second = labels_second;
end

if any(strcmpi(second_pass_mode, {'auto', 'target_linear_block', 'target_linear_pairs', 'target_linear_spatial_cvx'}))
    [xobs_fit_phys, ~] = build_mrdmd_xobs_for_labels( ...
        list_w, list_b, list_t_start, list_bin_widths, ...
        dt, m, mode_labels, fit_start_idx, fit_end_idx, list_anchor_idx);

    y_true_test_phys = xobs_test(:, target_cols) .* sigma_xobs(target_cols) + mu_xobs(target_cols);
    t_fit_feature = normalized_absolute_time_feature(fit_start_idx, fit_end_idx, fit_start_idx, fit_end_idx);
    t_test_feature = normalized_absolute_time_feature(test_start_idx, test_end_idx, fit_start_idx, fit_end_idx);

    if strcmpi(second_pass_mode, 'target_linear_spatial_cvx')
        [K_spatial_block, spatial_cvx_info] = fit_target_linear_spatial_cvx_model( ...
            xobs_fit_phys(:, target_cols), dt, list_modes, mode_labels, target_cols, ...
            target_linear_ridge_alpha, target_linear_include_constant, target_linear_include_time, ...
            t_fit_feature, spatial_cvx_beta, spatial_cvx_lambda_l1, spatial_cvx_sample_stride, ...
            spatial_cvx_amp_beta, X_fit_background, X_raw_fit);

        y_spatial_block_phys = predict_target_linear_one_step(y_true_test_phys, K_spatial_block, dt, ...
            target_linear_include_constant, target_linear_include_time, t_test_feature);

        [X_spatial_block, err_spatial_block, corr_spatial_block] = score_l3_candidate( ...
            y_spatial_block_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X);

        candidate_report(end+1) = make_candidate_report('target_linear_spatial_cvx', corr_spatial_block, err_spatial_block);

        best_name = 'target_linear_spatial_cvx';
        best_y_recreate_phys = y_spatial_block_phys;
        best_X_l3_wsindy = X_spatial_block;
        best_l3_recreate_error = err_spatial_block;
        best_l3_corr = corr_spatial_block;
        [best_w_second, best_labels_second] = make_linear_equations_for_targets( ...
            K_spatial_block, target_cols, mode_labels, 1:length(target_cols), w_second, labels_second, ...
            target_linear_include_constant, target_linear_include_time);
    end

    if any(strcmpi(second_pass_mode, {'auto', 'target_linear_block'}))
        y_block_phys = y_recreate_phys;
        K_block = fit_target_linear_model(xobs_fit_phys(:, target_cols), dt, ...
            target_linear_ridge_alpha, target_linear_include_constant, target_linear_include_time, t_fit_feature);
        y_block_phys(:, :) = predict_target_linear_one_step(y_true_test_phys, K_block, dt, ...
            target_linear_include_constant, target_linear_include_time, t_test_feature);
        [y_block_phys(:, :), calibration_report_block] = calibrate_target_linear_prediction( ...
            xobs_fit_phys(:, target_cols), y_block_phys(:, :), K_block, dt, ...
            target_linear_include_constant, target_linear_include_time, ...
            t_fit_feature, target_linear_calibrate_amplitude, mode_labels(target_cols));

        [X_block, err_block, corr_block] = score_l3_candidate( ...
            y_block_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X);

        candidate_report(end+1) = make_candidate_report('target_linear_block', corr_block, err_block);

        if is_better_l3_candidate(corr_block, err_block, best_l3_corr, best_l3_recreate_error)
            best_name = 'target_linear_block';
            best_y_recreate_phys = y_block_phys;
            best_X_l3_wsindy = X_block;
            best_l3_recreate_error = err_block;
            best_l3_corr = corr_block;
            [best_w_second, best_labels_second] = make_linear_equations_for_targets( ...
                K_block, target_cols, mode_labels, 1:length(target_cols), w_second, labels_second, ...
                target_linear_include_constant, target_linear_include_time);
            best_calibration_report = calibration_report_block;
        end
    end

    if any(strcmpi(second_pass_mode, {'auto', 'target_linear_pairs'})) && length(target_cols) >= 2
        for aa = 1:(length(target_cols)-1)
            for bb = (aa+1):length(target_cols)
                pair_idx = [aa bb];
                y_pair_phys = y_recreate_phys;
                K_pair = fit_target_linear_model(xobs_fit_phys(:, target_cols(pair_idx)), dt, ...
                    target_linear_ridge_alpha, target_linear_include_constant, target_linear_include_time, t_fit_feature);
                y_pair_phys(:, pair_idx) = predict_target_linear_one_step(y_true_test_phys(:, pair_idx), K_pair, dt, ...
                    target_linear_include_constant, target_linear_include_time, t_test_feature);
                [y_pair_phys(:, pair_idx), calibration_report_pair] = calibrate_target_linear_prediction( ...
                    xobs_fit_phys(:, target_cols(pair_idx)), y_pair_phys(:, pair_idx), K_pair, dt, ...
                    target_linear_include_constant, target_linear_include_time, ...
                    t_fit_feature, target_linear_calibrate_amplitude, mode_labels(target_cols(pair_idx)));

                [X_pair, err_pair, corr_pair] = score_l3_candidate( ...
                    y_pair_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X);

                pair_name = make_pair_candidate_name(mode_labels, target_cols(pair_idx));
                candidate_report(end+1) = make_candidate_report(pair_name, corr_pair, err_pair);

                if is_better_l3_candidate(corr_pair, err_pair, best_l3_corr, best_l3_recreate_error)
                    best_name = pair_name;
                    best_y_recreate_phys = y_pair_phys;
                    best_X_l3_wsindy = X_pair;
                    best_l3_recreate_error = err_pair;
                    best_l3_corr = corr_pair;
                    [best_w_second, best_labels_second] = make_linear_equations_for_targets( ...
                        K_pair, target_cols, mode_labels, pair_idx, w_second, labels_second, ...
                        target_linear_include_constant, target_linear_include_time);
                    best_calibration_report = calibration_report_pair;
                end
            end
        end
    end
end

y_recreate_phys = best_y_recreate_phys;
X_l3_wsindy = best_X_l3_wsindy;
l3_recreate_error = best_l3_recreate_error;
l3_corr = best_l3_corr;
w_second = best_w_second;
labels_second = best_labels_second;

if ~strcmpi(best_name, 'wsindy')
    y_recreate = (y_recreate_phys - mu_xobs(target_cols)) ./ sigma_xobs(target_cols);
end

X_raw_true = X(:, test_start_idx:test_end_idx);

% Use the official mrDMD reconstruction as the full modal truth. This
% includes every configured mrDMD level, not just the WSINDy input/target
% levels. The WSINDy model still only replaces the selected target-level
% modes.
X_inside_true = real(X_mrdmd_full(:, test_start_idx:test_end_idx));
X_inside_pred = real(X_inside_true - X_l3_true + X_l3_wsindy);

inside_full_error = zeros(1, recreate_steps);
raw_l3_error = zeros(1, recreate_steps);
raw_full_error = zeros(1, recreate_steps);
raw_full_frame_error = zeros(1, recreate_steps);
raw_full_corr = zeros(1, recreate_steps);
raw_full_cos = zeros(1, recreate_steps);
raw_frame_norm = zeros(1, recreate_steps);
pred_frame_norm = zeros(1, recreate_steps);
mrdmd_full_corr = zeros(1, recreate_steps);
mrdmd_raw_error = zeros(1, recreate_steps);

for k = 1:recreate_steps
    raw_frame_norm(k) = norm(X_raw_true(:,k));
    pred_frame_norm(k) = norm(X_inside_pred(:,k));
    inside_full_error(k) = norm(X_inside_true(:,k) - X_inside_pred(:,k)) / ...
                            (global_scale + eps);
    raw_l3_error(k) = norm(X_raw_true(:,k) - X_l3_wsindy(:,k)) / ...
                            (global_scale + eps);
    raw_full_error(k) = norm(X_raw_true(:,k) - X_inside_pred(:,k)) / ...
                            (global_scale + eps);
    raw_full_frame_error(k) = norm(X_raw_true(:,k) - X_inside_pred(:,k)) / ...
                            (raw_frame_norm(k) + eps);
    raw_full_corr(k) = centered_spatial_corr(X_raw_true(:,k), X_inside_pred(:,k));
    raw_full_cos(k) = real(dot(X_raw_true(:,k), X_inside_pred(:,k))) / ...
                            ((raw_frame_norm(k) + eps) * (pred_frame_norm(k) + eps));
    mrdmd_full_corr(k) = centered_spatial_corr(X_inside_true(:,k), X_inside_pred(:,k));
    mrdmd_raw_error(k) = norm(X_raw_true(:,k) - X_inside_true(:,k)) / ...
                            (global_scale + eps);
end

mean_l3_err = mean(l3_recreate_error) * 100;
max_l3_err = max(l3_recreate_error) * 100;
mean_full_err = mean(raw_full_frame_error) * 100;
max_full_err = max(raw_full_frame_error) * 100;

if nargout > 4
    details.w_second = w_second;
    details.labels_second = labels_second;
    details.mode_labels = mode_labels;
    details.target_cols = target_cols;
    details.mu_xobs = mu_xobs;
    details.sigma_xobs = sigma_xobs;
    details.xobs_test = xobs_test;
    details.tobs_test = tobs_test;
    details.test_indices = test_start_idx:test_end_idx;
    details.y_recreate = y_recreate;
    details.y_recreate_phys = y_recreate_phys;
    details.X_l3_true = X_l3_true;
    details.X_l3_wsindy = X_l3_wsindy;
    details.X_inside_true = X_inside_true;
    details.X_inside_pred = X_inside_pred;
    details.X_raw_true = X_raw_true;
    details.l3_recreate_error = l3_recreate_error;
    details.l3_spatial_corr = l3_corr;
    details.raw_l3_error = raw_l3_error;
    details.inside_full_error = inside_full_error;
    details.raw_full_error = raw_full_error;
    details.raw_full_frame_error = raw_full_frame_error;
    details.raw_full_corr = raw_full_corr;
    details.raw_full_cos = raw_full_cos;
    details.raw_frame_norm = raw_frame_norm;
    details.pred_frame_norm = pred_frame_norm;
    details.mrdmd_full_corr = mrdmd_full_corr;
    details.mrdmd_raw_error = mrdmd_raw_error;
    details.second_pass_mode = second_pass_mode;
    details.selected_second_pass_model = best_name;
    details.candidate_report = candidate_report;
    details.target_linear_ridge_alpha = target_linear_ridge_alpha;
    details.target_linear_include_constant = target_linear_include_constant;
    details.target_linear_include_time = target_linear_include_time;
    details.target_linear_calibrate_amplitude = target_linear_calibrate_amplitude;
    details.target_linear_calibration_report = best_calibration_report;
    details.spatial_cvx_beta = spatial_cvx_beta;
    details.spatial_cvx_amp_beta = spatial_cvx_amp_beta;
    details.spatial_cvx_lambda_l1 = spatial_cvx_lambda_l1;
    details.spatial_cvx_sample_stride = spatial_cvx_sample_stride;
    details.spatial_cvx_info = spatial_cvx_info;
end

end

function y_pred = predict_wsindy_one_step(xobs, tobs, target_cols, ...
    w_second, labels_second, mode_labels, dt)
%PREDICT_WSINDY_ONE_STEP Teacher-forced one-step target prediction.

y_pred = xobs(:, target_cols);

for q = 1:(size(xobs, 1)-1)
    y_now = xobs(q, target_cols).';
    dydt_now = mrdmd_wsindy_inside_rhs( ...
        tobs(q), y_now, ...
        w_second, labels_second, mode_labels, target_cols, ...
        xobs, tobs);

    y_pred(q+1, :) = xobs(q, target_cols) + dt * dydt_now.';
end

end

function [X_l3_candidate, l3_error, l3_corr] = score_l3_candidate( ...
    y_candidate_phys, list_modes, mode_labels, target_cols, X_l3_true, global_scale, n, X)
%SCORE_L3_CANDIDATE Rebuild selected L3 field and score shape plus magnitude.

recreate_steps = size(y_candidate_phys, 1);
X_l3_candidate = build_l3_contribution_from_amplitudes( ...
    y_candidate_phys, list_modes, mode_labels, target_cols, recreate_steps, n, X);

l3_error = zeros(1, recreate_steps);
l3_corr = zeros(1, recreate_steps);

for k = 1:recreate_steps
    l3_error(k) = norm(X_l3_true(:,k) - X_l3_candidate(:,k)) / (global_scale + eps);
    l3_corr(k) = centered_spatial_corr(X_l3_true(:,k), X_l3_candidate(:,k));
end

end

function X_l3 = build_l3_contribution_from_amplitudes( ...
    y_phys, list_modes, mode_labels, target_cols, recreate_steps, n, X)
%BUILD_L3_CONTRIBUTION_FROM_AMPLITUDES Convert target amplitudes to fields.

X_l3 = zeros(n, recreate_steps, 'like', X);

for kk = 1:length(target_cols)
    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);
    phi = list_modes{lev,bin}(:,mode_idx);

    for q = 1:recreate_steps
        X_l3(:,q) = X_l3(:,q) + phi * y_phys(q,kk);
    end
end

X_l3 = real(X_l3);

end

function K = fit_target_linear_model(A_fit, dt, alpha, include_constant, include_time, t_feature)
%FIT_TARGET_LINEAR_MODEL Dense target-only linear derivative model.

if size(A_fit, 1) < 2
    K = zeros(size(A_fit, 2) + include_constant + include_time, size(A_fit, 2));
    return;
end

dA = (A_fit(2:end,:) - A_fit(1:end-1,:)) / dt;
Theta = build_target_linear_features(A_fit(1:end-1,:), include_constant, include_time, t_feature(1:end-1));

if alpha > 0
    K = (Theta.'*Theta + alpha*eye(size(Theta,2))) \ (Theta.'*dA);
else
    K = Theta \ dA;
end

end

function [K, info] = fit_target_linear_spatial_cvx_model( ...
    A_fit, dt, list_modes, mode_labels, target_cols, alpha_ridge, ...
    include_constant, include_time, t_feature, beta_spatial, lambda_l1, sample_stride, ...
    beta_amp, X_background_fit, X_raw_fit)
%FIT_TARGET_LINEAR_SPATIAL_CVX_MODEL Dense target block with full-frame one-step CVX loss.

if nargin < 7 || isempty(alpha_ridge)
    alpha_ridge = 0;
end
if nargin < 8 || isempty(include_constant)
    include_constant = false;
end
if nargin < 9 || isempty(include_time)
    include_time = false;
end
if nargin < 10 || isempty(t_feature)
    t_feature = zeros(size(A_fit, 1), 1);
end
if nargin < 11 || isempty(beta_spatial)
    beta_spatial = 1;
end
if nargin < 12 || isempty(lambda_l1)
    lambda_l1 = 0;
end
if nargin < 13 || isempty(sample_stride)
    sample_stride = 1;
end
if nargin < 14 || isempty(beta_amp)
    beta_amp = 1;
end
if nargin < 15
    X_background_fit = [];
end
if nargin < 16
    X_raw_fit = [];
end
A_fit = double(real(A_fit));
t_feature = double(real(t_feature(:)));

if exist('cvx_begin', 'file') ~= 2
    error(['CVX is not on the MATLAB path. Install/add CVX first, then rerun ', ...
        'target_linear_spatial_cvx.']);
end

if size(A_fit, 1) < 2
    K = zeros(size(A_fit, 2) + include_constant + include_time, size(A_fit, 2));
    info.status = 'not_enough_samples';
    info.optval = NaN;
    return;
end

sample_idx = 1:sample_stride:(size(A_fit, 1)-1);
Theta_all = build_target_linear_features(A_fit(1:end-1,:), include_constant, include_time, t_feature(1:end-1));
dA_all = (A_fit(2:end,:) - A_fit(1:end-1,:)) / dt;

Theta_sample = build_target_linear_features(A_fit(sample_idx,:), include_constant, include_time, t_feature(sample_idx));
A_now_sample = A_fit(sample_idx, :);
A_next_true_sample = A_fit(sample_idx+1, :);

Phi = double(build_target_phi_real(list_modes, mode_labels, target_cols));
X_next_true = A_next_true_sample * Phi.';
use_full_frame_spatial = ~isempty(X_background_fit) && ~isempty(X_raw_fit);

if use_full_frame_spatial
    X_background_fit = double(real(X_background_fit));
    X_raw_fit = double(real(X_raw_fit));
    X_next_background = X_background_fit(:, sample_idx+1).';
    X_next_raw = X_raw_fit(:, sample_idx+1).';
else
    X_next_background = zeros(size(X_next_true));
    X_next_raw = X_next_true;
end

frame_energy = vecnorm(X_next_raw, 2, 2).^2;
if max(frame_energy) > 0
    frame_weights = frame_energy / mean(frame_energy(frame_energy > 0));
else
    frame_weights = ones(size(frame_energy));
end
sqrt_weights = sqrt(frame_weights);

spatial_scale = max(norm(X_next_raw, 'fro')^2, eps);
amp_scale = max(norm(A_next_true_sample, 'fro')^2, eps);
derivative_scale = max(norm(dA_all, 'fro')^2, eps);
num_features = size(Theta_all, 2);
num_targets = size(A_fit, 2);
num_deriv_samples = size(dA_all, 1);
num_samples = length(sample_idx);
num_spatial_states = size(Phi, 1);

cvx_begin quiet
    variable K_cvx(num_features, num_targets)

    expression dA_pred(num_deriv_samples, num_targets)
    expression A_next_pred(num_samples, num_targets)
    expression X_next_pred(num_samples, num_spatial_states)
    expression X_next_full_pred(num_samples, num_spatial_states)
    expression spatial_residual(num_samples, num_spatial_states)
    expression amp_residual(num_samples, num_targets)

    dA_pred = Theta_all * K_cvx;
    A_next_pred = A_now_sample + dt * (Theta_sample * K_cvx);
    X_next_pred = A_next_pred * Phi.';
    X_next_full_pred = X_next_background + X_next_pred;
    spatial_residual = (X_next_full_pred - X_next_raw) .* repmat(sqrt_weights, 1, num_spatial_states);
    amp_residual = (A_next_pred - A_next_true_sample) .* repmat(sqrt_weights, 1, num_targets);

    minimize( ...
        sum_square(dA_pred(:) - dA_all(:)) / derivative_scale + ...
        beta_amp * sum_square(amp_residual(:)) / amp_scale + ...
        beta_spatial * sum_square(spatial_residual(:)) / spatial_scale + ...
        alpha_ridge * sum_square(K_cvx(:)) + ...
        lambda_l1 * norm(K_cvx(:), 1))
cvx_end

K = K_cvx;
info.status = cvx_status;
info.optval = cvx_optval;
info.beta_spatial = beta_spatial;
info.beta_amp = beta_amp;
info.lambda_l1 = lambda_l1;
info.alpha_ridge = alpha_ridge;
info.sample_stride = sample_stride;
info.sample_idx = sample_idx;
info.mean_frame_weight = mean(frame_weights);
info.max_frame_weight = max(frame_weights);
info.use_full_frame_spatial = use_full_frame_spatial;
A_next_pred_num = A_now_sample + dt * (Theta_sample * K);
X_next_pred_num = A_next_pred_num * Phi.';
X_next_full_pred_num = X_next_background + X_next_pred_num;
info.derivative_loss = norm(Theta_all * K - dA_all, 'fro')^2 / derivative_scale;
info.amplitude_loss = norm((A_next_pred_num - A_next_true_sample) .* ...
    repmat(sqrt_weights, 1, num_targets), 'fro')^2 / amp_scale;
info.spatial_loss = norm((X_next_full_pred_num - X_next_raw) .* ...
    repmat(sqrt_weights, 1, num_spatial_states), 'fro')^2 / spatial_scale;

end

function Phi = build_target_phi_real(list_modes, mode_labels, target_cols)
%BUILD_TARGET_PHI_REAL Real selected-target spatial basis.

first_label = mode_labels{target_cols(1)};
[lev, bin, ~] = parse_mode_label(first_label);
n = size(list_modes{lev, bin}, 1);
Phi = zeros(n, length(target_cols));

for kk = 1:length(target_cols)
    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);
    Phi(:, kk) = real(list_modes{lev, bin}(:, mode_idx));
end

end

function A_pred = predict_target_linear_one_step(A_true, K, dt, include_constant, include_time, t_feature)
%PREDICT_TARGET_LINEAR_ONE_STEP Teacher-forced one-step modal prediction.

A_pred = A_true;

for k = 1:(size(A_true, 1)-1)
    theta_now = build_target_linear_features(A_true(k,:), include_constant, include_time, t_feature(k));
    A_pred(k+1,:) = A_true(k,:) + dt * (theta_now * K);
end

end

function [A_cal, report] = calibrate_target_linear_prediction( ...
    A_fit_true, A_pred_test, K, dt, include_constant, include_time, t_fit_feature, do_calibrate, labels)
%CALIBRATE_TARGET_LINEAR_PREDICTION Scale/bias target-linear outputs when shape is already right.

amp_corr_threshold = 0.95;
rel_err_threshold = 0.15;

A_cal = A_pred_test;
A_fit_pred = predict_target_linear_one_step(A_fit_true, K, dt, include_constant, include_time, t_fit_feature);
report = struct('label', {}, 'applied', {}, 'fit_amp_corr', {}, 'fit_rel_err', {}, 'scale', {}, 'bias', {});

for jj = 1:size(A_fit_true, 2)
    true_fit = A_fit_true(:, jj);
    pred_fit = A_fit_pred(:, jj);

    fit_corr = amplitude_corr(true_fit, pred_fit);
    fit_rel_err = norm(true_fit - pred_fit) / (norm(true_fit) + eps);

    scale = 1;
    bias = 0;
    applied = false;

    if do_calibrate && fit_corr > amp_corr_threshold && fit_rel_err > rel_err_threshold
        pred_centered = pred_fit - mean(pred_fit);
        true_centered = true_fit - mean(true_fit);
        scale = (pred_centered(:).' * true_centered(:)) / ...
            (pred_centered(:).' * pred_centered(:) + eps);
        bias = mean(true_fit) - scale * mean(pred_fit);
        A_cal(:, jj) = scale * A_pred_test(:, jj) + bias;
        applied = true;
    end

    report(end+1).label = labels{jj};
    report(end).applied = applied;
    report(end).fit_amp_corr = fit_corr;
    report(end).fit_rel_err = fit_rel_err;
    report(end).scale = scale;
    report(end).bias = bias;
end

end

function rho = amplitude_corr(a, b)
%AMPLITUDE_CORR Correlation for 1D modal amplitudes with flat-signal guard.

a = real(a(:));
b = real(b(:));
a = a - mean(a);
b = b - mean(b);

denom = norm(a) * norm(b);
if denom <= eps
    rho = 0;
else
    rho = dot(a, b) / denom;
end

end

function Theta = build_target_linear_features(A, include_constant, include_time, t_feature)
%BUILD_TARGET_LINEAR_FEATURES Restricted target-linear feature matrix.

Theta = A;

if include_constant
    Theta = [ones(size(A, 1), 1), Theta];
end

if include_time
    Theta = [Theta, t_feature(:)];
end

end

function t_feature = normalized_absolute_time_feature(start_idx, end_idx, fit_start_idx, fit_end_idx)
%NORMALIZED_ABSOLUTE_TIME_FEATURE Time relative to the fit window.

frames = (start_idx:end_idx).';
denom = max(fit_end_idx - fit_start_idx, 1);
t_feature = 2 * (frames - fit_start_idx) / denom - 1;

end

function tf = is_better_l3_candidate(corr_new, err_new, corr_best, err_best)
%IS_BETTER_L3_CANDIDATE Prefer selected-L3 shape, then selected-L3 error.

mean_corr_new = mean(corr_new);
mean_corr_best = mean(corr_best);
mean_err_new = mean(err_new);
mean_err_best = mean(err_best);

tf = mean_corr_new > mean_corr_best + 1e-6 || ...
    (abs(mean_corr_new - mean_corr_best) <= 1e-6 && mean_err_new < mean_err_best);

end

function report = make_candidate_report(name, l3_corr, l3_error)
%MAKE_CANDIDATE_REPORT Compact diagnostics for model selection.

report.name = name;
report.mean_l3_corr = mean(l3_corr);
report.mean_l3_error = mean(l3_error) * 100;
report.max_l3_error = max(l3_error) * 100;

end

function name = make_pair_candidate_name(mode_labels, pair_cols)
%MAKE_PAIR_CANDIDATE_NAME Human-readable pair model name.

label_a = regexprep(mode_labels{pair_cols(1)}, '\s+', '_');
label_b = regexprep(mode_labels{pair_cols(2)}, '\s+', '_');
name = sprintf('target_linear_pair_%s__%s', label_a, label_b);

end

function [w_out, labels_out] = make_linear_equations_for_targets( ...
    K, target_cols, mode_labels, local_target_idx, w_base, labels_base, include_constant, include_time)
%MAKE_LINEAR_EQUATIONS_FOR_TARGETS Replace selected equations with target-linear equations.

w_out = w_base;
labels_out = labels_base;
linear_labels = mode_labels(target_cols(local_target_idx));
if include_constant
    linear_labels = [{'1'}, linear_labels];
end
if include_time
    linear_labels = [linear_labels, {'t'}];
end

for jj = 1:length(local_target_idx)
    eq_idx = local_target_idx(jj);
    w_out{eq_idx} = K(:, jj);
    labels_out{eq_idx} = linear_labels;
end

end

function X_l3_true = build_true_l3_contribution_from_modes( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    list_anchor_idx, mode_labels, target_cols, test_start_idx, recreate_steps, n, X)
%BUILD_TRUE_L3_CONTRIBUTION_FROM_MODES Match reconstruct_mrdmd for selected modes.

X_l3_true = zeros(n, recreate_steps, 'like', X);

for kk = 1:length(target_cols)
    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);

    modes = list_modes{lev, bin};
    eigs_slow = list_w{lev, bin};
    b = list_b{lev, bin};
    t_start = list_t_start(lev, bin);
    bin_width = list_bin_widths(lev, bin);
    t_end = t_start + bin_width - 1;
    anchor_idx = list_anchor_idx(lev, bin);
    anchor_frame = t_start + anchor_idx - 1;

    for q = 1:recreate_steps
        absolute_frame = test_start_idx + q - 1;

        if absolute_frame < t_start || absolute_frame > t_end || absolute_frame < anchor_frame
            continue;
        end

        rel_time = absolute_frame - anchor_frame;
        amp = b(mode_idx) * eigs_slow(mode_idx)^rel_time;
        X_l3_true(:,q) = X_l3_true(:,q) + modes(:,mode_idx) * amp;
    end
end

X_l3_true = real(X_l3_true);

end

function rho = centered_spatial_corr(a, b)
%CENTERED_SPATIAL_CORR Shape agreement, ignoring constant offset.

a = real(a(:));
b = real(b(:));
a = a - mean(a);
b = b - mean(b);
rho = dot(a, b) / ((norm(a) * norm(b)) + eps);

end
