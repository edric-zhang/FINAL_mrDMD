function [mean_l3_err, max_l3_err, mean_full_err, max_full_err, details] = evaluate_wsindy_parameter_set( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    X, dt, m, n, input_levels, target_level, ...
    fit_start_idx, fit_end_idx, test_start_idx, test_end_idx, ...
    top_input_modes_per_level, top_target_modes, ...
    lambda1, lambda2, gamma, max_terms_per_equation, max_quadratic_base_terms)
%EVALUATE_WSINDY_PARAMETER_SET Train and score one WSINDy parameter set.

[w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs] = run_mrdmd_wsindy( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    dt, m, input_levels, target_level, fit_start_idx, fit_end_idx, ...
    top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
    max_terms_per_equation, max_quadratic_base_terms);

test_steps = test_end_idx - test_start_idx + 1;
recreate_steps = test_steps;

[xobs_test, tobs_test] = build_mrdmd_xobs_for_labels( ...
    list_w, list_b, list_t_start, list_bin_widths, ...
    dt, m, mode_labels, test_start_idx, test_end_idx);

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

X_l3_true = zeros(n, recreate_steps, 'like', X);
X_l3_wsindy = zeros(n, recreate_steps, 'like', X);

% Build true L3 contribution from actual modal amplitudes in xobs
for kk = 1:length(target_cols)

    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);

    phi = list_modes{lev,bin}(:,mode_idx);

    for q = 1:recreate_steps
        true_amp_phys = xobs_test(q,target_cols(kk)) .* sigma_xobs(target_cols(kk)) + mu_xobs(target_cols(kk));
        X_l3_true(:,q) = X_l3_true(:,q) + phi * true_amp_phys;
        X_l3_wsindy(:,q) = X_l3_wsindy(:,q) + phi * y_recreate_phys(q,kk);
    end
end

X_l3_true = real(X_l3_true);
X_l3_wsindy = real(X_l3_wsindy);

l3_recreate_error = zeros(1, recreate_steps);

for k = 1:recreate_steps
    l3_recreate_error(k) = norm(X_l3_true(:,k) - X_l3_wsindy(:,k)) / ...
                            (norm(X_l3_true(:,k)) + eps);
end

% FULL INSIDE-BIN CFD RECREATION
% Add the unchanged L1/L2 ancestor contributions inside the same bin,
% then add the WSINDy-recreated L3 contribution.
X_inside_pred = zeros(n, recreate_steps, 'like', X);

% Add L1/L2 true DMD contributions over the held-out test frames
for lev = 1:2
    J = 2^(lev-1);

    for bin = 1:J
        if isempty(list_modes{lev,bin})
            continue;
        end

        t_start = list_t_start(lev,bin);
        bin_width = list_bin_widths(lev,bin);
        t_end = t_start + bin_width - 1;

        if t_start > test_end_idx || t_end < test_start_idx
            continue;
        end

        modes = list_modes{lev,bin};
        eigs_slow = list_w{lev,bin};
        b = list_b{lev,bin};

        for q = 1:recreate_steps
            absolute_frame = test_start_idx + q - 1;

            if absolute_frame < t_start || absolute_frame > t_end
                continue;
            end

            rel_time = absolute_frame - t_start;
            time_powers = eigs_slow .^ rel_time;

            X_inside_pred(:,q) = X_inside_pred(:,q) + modes * (b .* time_powers);
        end
    end
end

% Add WSINDy-recreated L3 contribution
X_inside_pred = real(X_inside_pred + X_l3_wsindy);
X_inside_true = X(:, test_start_idx:test_end_idx);

inside_full_error = zeros(1, recreate_steps);

for k = 1:recreate_steps
    inside_full_error(k) = norm(X_inside_true(:,k) - X_inside_pred(:,k)) / ...
                            (norm(X_inside_true(:,k)) + eps);
end

mean_l3_err = mean(l3_recreate_error) * 100;
max_l3_err = max(l3_recreate_error) * 100;
mean_full_err = mean(inside_full_error) * 100;
max_full_err = max(inside_full_error) * 100;

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
    details.l3_recreate_error = l3_recreate_error;
    details.inside_full_error = inside_full_error;
end

end
