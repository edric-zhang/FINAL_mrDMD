%RUN_SPATIAL_CVX_REFIT_FROM_DETAILS Experimental selected-L3 spatial CVX refit.
%
% Run this after main/run_single_experiment.m, once CVX is installed and on
% the MATLAB path. It refits the currently selected equation supports using
% derivative loss plus one-step selected-L3 spatial loss on the fit window.

if ~exist('details', 'var') || ~exist('mrdmd', 'var') || ~exist('data', 'var') || ~exist('cfg', 'var')
    error('Run main/run_single_experiment.m first so details, mrdmd, data, and cfg exist.');
end

beta_spatial = 1;
lambda_l1 = 0;
sample_stride = 5;

[xobs_fit_raw, tobs_fit] = build_mrdmd_xobs_for_labels( ...
    mrdmd.list_w, mrdmd.list_b, mrdmd.list_t_start, mrdmd.list_bin_widths, ...
    data.dt, data.m, details.mode_labels, ...
    cfg.frames.fit_start_idx, cfg.frames.fit_end_idx, mrdmd.list_anchor_idx);

xobs_fit = (xobs_fit_raw - details.mu_xobs) ./ details.sigma_xobs;

[w_spatial, spatial_cvx_info] = run_mrdmd_wsindy_spatial_cvx( ...
    details.labels_second, details.mode_labels, details.target_cols, ...
    xobs_fit, tobs_fit, details.mu_xobs, details.sigma_xobs, ...
    mrdmd.list_modes, beta_spatial, lambda_l1, sample_stride);

fprintf('\n================ SPATIAL-CVX REFIT EQUATIONS ================\n');
print_wsindy_equations(w_spatial, details.labels_second, details.mode_labels, details.target_cols);

disp(spatial_cvx_info);
