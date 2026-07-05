function [w_spatial, info] = run_mrdmd_wsindy_spatial_cvx( ...
    labels_second, mode_labels, target_cols, xobs, tobs, ...
    mu_xobs, sigma_xobs, list_modes, ...
    beta_spatial, lambda_l1, sample_stride)
%RUN_MRDMD_WSINDY_SPATIAL_CVX Refit WSINDy coefficients with one-step L3 spatial loss.
%
% This is an experimental CVX-based refinement. It keeps the existing
% WSINDy dictionary/support in labels_second, but refits coefficients using:
%
%   derivative loss
% + beta_spatial * selected-L3 one-step spatial loss
% + lambda_l1 * coefficient L1 penalty
%
% The one-step spatial term uses:
%   a_pred(k+1) = a(k) + dt * f(a(k), inputs(k))
%   X_L3_pred(k+1) = Phi_target * a_pred(k+1)
%
% This function requires CVX on the MATLAB path.

if nargin < 10 || isempty(beta_spatial)
    beta_spatial = 1;
end
if nargin < 11 || isempty(lambda_l1)
    lambda_l1 = 0;
end
if nargin < 12 || isempty(sample_stride)
    sample_stride = 5;
end

if exist('cvx_begin', 'file') ~= 2
    error(['CVX is not on the MATLAB path. Install/add CVX first, then rerun ', ...
        'run_mrdmd_wsindy_spatial_cvx.']);
end

dt = tobs(2) - tobs(1);
num_targets = length(target_cols);
sample_idx = 1:sample_stride:(size(xobs, 1)-1);

Theta_cells = cell(num_targets, 1);
fd_targets = cell(num_targets, 1);
coef_offsets = zeros(num_targets, 1);
coef_lengths = zeros(num_targets, 1);

total_coef = 0;
for kk = 1:num_targets
    labels = labels_second{kk};
    Theta_cells{kk} = build_theta_for_labels(labels, sample_idx, xobs, tobs, ...
        mode_labels, target_cols);
    fd_targets{kk} = (xobs(sample_idx+1, target_cols(kk)) - ...
        xobs(sample_idx, target_cols(kk))) / dt;

    coef_offsets(kk) = total_coef;
    coef_lengths(kk) = length(labels);
    total_coef = total_coef + coef_lengths(kk);
end

Phi = build_target_phi(list_modes, mode_labels, target_cols);
A_now = xobs(sample_idx, target_cols);
A_next_true = xobs(sample_idx+1, target_cols);

A_now_phys = A_now .* sigma_xobs(target_cols) + mu_xobs(target_cols);
A_next_true_phys = A_next_true .* sigma_xobs(target_cols) + mu_xobs(target_cols);
X_next_true = A_next_true_phys * Phi.';

global_spatial_scale = max(norm(X_next_true, 'fro'), eps);
num_samples = length(sample_idx);

cvx_begin quiet
    variable coef_all(total_coef)

    expression deriv_loss
    expression spatial_loss
    expression A_next_pred_norm(num_samples, num_targets)
    expression A_next_pred_phys(num_samples, num_targets)
    expression X_next_pred(num_samples, size(Phi, 1))

    deriv_loss = 0;

    for kk = 1:num_targets
        idx = (coef_offsets(kk)+1):(coef_offsets(kk)+coef_lengths(kk));
        wk = coef_all(idx);
        dydt_pred = Theta_cells{kk} * wk;

        deriv_loss = deriv_loss + sum_square(dydt_pred - fd_targets{kk});
        A_next_pred_norm(:, kk) = A_now(:, kk) + dt * dydt_pred;
    end

    for kk = 1:num_targets
        A_next_pred_phys(:, kk) = A_next_pred_norm(:, kk) * ...
            sigma_xobs(target_cols(kk)) + mu_xobs(target_cols(kk));
    end

    X_next_pred = A_next_pred_phys * Phi.';
    spatial_loss = sum_square(X_next_pred(:) - X_next_true(:)) / global_spatial_scale^2;

    minimize(deriv_loss + beta_spatial * spatial_loss + lambda_l1 * norm(coef_all, 1))
cvx_end

w_spatial = cell(num_targets, 1);
for kk = 1:num_targets
    idx = (coef_offsets(kk)+1):(coef_offsets(kk)+coef_lengths(kk));
    w_spatial{kk} = coef_all(idx);
end

info.status = cvx_status;
info.optval = cvx_optval;
info.beta_spatial = beta_spatial;
info.lambda_l1 = lambda_l1;
info.sample_stride = sample_stride;
info.sample_idx = sample_idx;
info.derivative_loss = evaluate_derivative_loss(w_spatial, Theta_cells, fd_targets);
info.spatial_loss = evaluate_spatial_loss(w_spatial, Theta_cells, A_now, ...
    A_next_true_phys, Phi, dt, mu_xobs(target_cols), sigma_xobs(target_cols));

end

function Theta = build_theta_for_labels(labels, sample_idx, xobs, tobs, mode_labels, target_cols)
%BUILD_THETA_FOR_LABELS Numeric library matrix for sampled one-step times.

Theta = zeros(length(sample_idx), length(labels));

for rr = 1:length(sample_idx)
    q = sample_idx(rr);
    y_now = xobs(q, target_cols).';

    for cc = 1:length(labels)
        Theta(rr, cc) = evaluate_inside_label(labels{cc}, tobs(q), y_now, ...
            mode_labels, target_cols, xobs, tobs);
    end
end

end

function Phi = build_target_phi(list_modes, mode_labels, target_cols)
%BUILD_TARGET_PHI Real selected-target spatial basis.

first_label = mode_labels{target_cols(1)};
[lev, bin, mode_idx] = parse_mode_label(first_label);
n = size(list_modes{lev, bin}, 1);
Phi = zeros(n, length(target_cols));

for kk = 1:length(target_cols)
    label = mode_labels{target_cols(kk)};
    [lev, bin, mode_idx] = parse_mode_label(label);
    Phi(:, kk) = real(list_modes{lev, bin}(:, mode_idx));
end

end

function loss = evaluate_derivative_loss(w_cells, Theta_cells, fd_targets)
%EVALUATE_DERIVATIVE_LOSS Numeric derivative loss after CVX solve.

loss = 0;
for kk = 1:length(w_cells)
    residual = Theta_cells{kk} * w_cells{kk} - fd_targets{kk};
    loss = loss + sum(residual.^2);
end

end

function loss = evaluate_spatial_loss(w_cells, Theta_cells, A_now, A_next_true_phys, ...
    Phi, dt, mu_target, sigma_target)
%EVALUATE_SPATIAL_LOSS Numeric selected-L3 one-step spatial loss after CVX solve.

A_next_pred_norm = A_now;

for kk = 1:length(w_cells)
    A_next_pred_norm(:, kk) = A_now(:, kk) + dt * (Theta_cells{kk} * w_cells{kk});
end

A_next_pred_phys = A_next_pred_norm .* sigma_target + mu_target;
X_pred = A_next_pred_phys * Phi.';
X_true = A_next_true_phys * Phi.';

loss = norm(X_pred - X_true, 'fro')^2 / max(norm(X_true, 'fro')^2, eps);

end
