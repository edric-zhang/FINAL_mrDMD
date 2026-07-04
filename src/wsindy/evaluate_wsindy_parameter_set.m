function [mean_l3_err, max_l3_err, mean_full_err, max_full_err, details] = evaluate_wsindy_parameter_set( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    X, dt, m, n, input_levels, target_level, ...
    fit_start_idx, fit_end_idx, test_start_idx, test_end_idx, ...
    top_input_modes_per_level, top_target_modes, ...
    lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms, list_anchor_idx)
%EVALUATE_WSINDY_PARAMETER_SET Train and score one WSINDy parameter set.

[w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs] = run_mrdmd_wsindy( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    dt, m, input_levels, target_level, fit_start_idx, fit_end_idx, ...
    top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
    max_terms_per_equation, max_quadratic_base_terms, list_anchor_idx);

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

X_l3_wsindy = zeros(n, recreate_steps, 'like', X);

% Build true L3 contribution from actual modal amplitudes in xobs
for kk = 1:length(target_cols)

    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);

    phi = list_modes{lev,bin}(:,mode_idx);

    for q = 1:recreate_steps
        X_l3_wsindy(:,q) = X_l3_wsindy(:,q) + phi * y_recreate_phys(q,kk);
    end
end

X_l3_true = build_true_l3_contribution_from_modes( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    list_anchor_idx, mode_labels, target_cols, test_start_idx, recreate_steps, n, X);
X_l3_wsindy = real(X_l3_wsindy);

l3_recreate_error = zeros(1, recreate_steps);

snapshot_norms = vecnorm(X);
valid_norms = snapshot_norms(snapshot_norms > 0);
global_scale = median(valid_norms);

for k = 1:recreate_steps
    l3_recreate_error(k) = norm(X_l3_true(:,k) - X_l3_wsindy(:,k)) / ...
                            (global_scale + eps);
end

X_raw_true = X(:, test_start_idx:test_end_idx);

% Use the official mrDMD reconstruction as the full modal truth. This
% includes every configured mrDMD level, not just the WSINDy input/target
% levels. The WSINDy model still only replaces the selected target-level
% modes.
mrdmd_tmp.list_w = list_w;
mrdmd_tmp.list_b = list_b;
mrdmd_tmp.list_modes = list_modes;
mrdmd_tmp.list_t_start = list_t_start;
mrdmd_tmp.list_bin_widths = list_bin_widths;
mrdmd_tmp.list_anchor_idx = list_anchor_idx;
mrdmd_tmp.L = size(list_modes, 1);

[X_mrdmd_full, ~] = reconstruct_mrdmd(mrdmd_tmp, X, 1, size(X, 2));
X_inside_true = real(X_mrdmd_full(:, test_start_idx:test_end_idx));
X_inside_pred = real(X_inside_true - X_l3_true + X_l3_wsindy);

inside_full_error = zeros(1, recreate_steps);
raw_l3_error = zeros(1, recreate_steps);
raw_full_error = zeros(1, recreate_steps);
raw_full_corr = zeros(1, recreate_steps);
mrdmd_full_corr = zeros(1, recreate_steps);
mrdmd_raw_error = zeros(1, recreate_steps);

for k = 1:recreate_steps
    inside_full_error(k) = norm(X_inside_true(:,k) - X_inside_pred(:,k)) / ...
                            (global_scale + eps);
    raw_l3_error(k) = norm(X_raw_true(:,k) - X_l3_wsindy(:,k)) / ...
                            (global_scale + eps);
    raw_full_error(k) = norm(X_raw_true(:,k) - X_inside_pred(:,k)) / ...
                            (global_scale + eps);
    raw_full_corr(k) = centered_spatial_corr(X_raw_true(:,k), X_inside_pred(:,k));
    mrdmd_full_corr(k) = centered_spatial_corr(X_inside_true(:,k), X_inside_pred(:,k));
    mrdmd_raw_error(k) = norm(X_raw_true(:,k) - X_inside_true(:,k)) / ...
                            (global_scale + eps);
end

mean_l3_err = mean(l3_recreate_error) * 100;
max_l3_err = max(l3_recreate_error) * 100;
mean_full_err = mean(raw_full_error) * 100;
max_full_err = max(raw_full_error) * 100;

if nargout > 4
    details.w_second = w_second;
    details.labels_second = labels_second;
    details.mode_labels = mode_labels;
    details.target_cols = target_cols;
    details.mu_xobs = mu_xobs;
    details.sigma_xobs = sigma_xobs;
    details.xobs_test = xobs_test;
    details.tobs_test = tobs_test;
    details.y_recreate = y_recreate;
    details.y_recreate_phys = y_recreate_phys;
    details.X_l3_true = X_l3_true;
    details.X_l3_wsindy = X_l3_wsindy;
    details.X_inside_true = X_inside_true;
    details.X_inside_pred = X_inside_pred;
    details.X_raw_true = X_raw_true;
    details.l3_recreate_error = l3_recreate_error;
    details.raw_l3_error = raw_l3_error;
    details.inside_full_error = inside_full_error;
    details.raw_full_error = raw_full_error;
    details.raw_full_corr = raw_full_corr;
    details.mrdmd_full_corr = mrdmd_full_corr;
    details.mrdmd_raw_error = mrdmd_raw_error;
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
