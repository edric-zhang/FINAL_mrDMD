function plot_raw_mrdmd_wsindy_comparison(X_raw_true, X_mrdmd_true, X_inside_pred, ...
    data, test_start_idx, recreate_steps)
%PLOT_RAW_MRDMD_WSINDY_COMPARISON Animate raw, true mrDMD, and WSINDy recreation.

figure('WindowStyle', 'normal', ...
    'Name', 'Raw vs mrDMD vs WSINDy Reconstruction', ...
    'NumberTitle', 'off', ...
    'Position', [80 100 1650 500]);

shared_clim = robust_three_way_clim(X_raw_true, X_mrdmd_true, X_inside_pred, data);

for k = 1:recreate_steps
    absolute_frame = test_start_idx + k - 1;

    ax1 = subplot(1,3,1);
    plot_state_field(ax1, X_raw_true(:,k), data, shared_clim, ...
        sprintf('Raw %s, Frame %d', data.plot_field, absolute_frame));

    ax2 = subplot(1,3,2);
    plot_state_field(ax2, X_mrdmd_true(:,k), data, shared_clim, ...
        sprintf('Full mrDMD %s, Frame %d', data.plot_field, absolute_frame));

    ax3 = subplot(1,3,3);
    plot_state_field(ax3, X_inside_pred(:,k), data, shared_clim, ...
        sprintf('Full WSINDy recreation %s, Frame %d', data.plot_field, absolute_frame));

    drawnow;
end

end

function clim_vals = robust_three_way_clim(X_raw_true, X_mrdmd_true, X_inside_pred, data)

all_values = [];
state_sets = {X_raw_true, X_mrdmd_true, X_inside_pred};

for ss = 1:numel(state_sets)
    X_set = state_sets{ss};
    for k = 1:size(X_set, 2)
        vals = state_to_field(X_set(:,k), data);
        vals = real(vals(:));
        vals = vals(isfinite(vals));
        all_values = [all_values; vals]; %#ok<AGROW>
    end
end

if isempty(all_values)
    clim_vals = [-1 1];
    return;
end

if strcmpi(data.plot_field, 'speed')
    hi = prctile(abs(all_values), 99);
    if hi <= 0 || ~isfinite(hi), hi = 1; end
    clim_vals = [0 hi];
else
    bound = prctile(abs(all_values), 99);
    if bound <= 0 || ~isfinite(bound), bound = 1; end
    clim_vals = [-bound bound];
end

end
