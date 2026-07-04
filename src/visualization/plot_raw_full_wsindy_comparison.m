function plot_raw_full_wsindy_comparison(X_raw_true, X_inside_pred, data, ...
    test_start_idx, recreate_steps)
%PLOT_RAW_FULL_WSINDY_COMPARISON Animate raw data vs WSINDy full recreation.

figure('WindowStyle', 'normal', ...
    'Name', 'Raw vs Full WSINDy Reconstruction', ...
    'NumberTitle', 'off', ...
    'Position', [100 100 1200 500]);

raw_clim = robust_raw_clim(X_raw_true, data);

for k = 1:recreate_steps
    absolute_frame = test_start_idx + k - 1;

    ax1 = subplot(1,2,1);
    plot_state_field(ax1, X_raw_true(:,k), data, raw_clim, ...
        sprintf('Raw %s, Frame %d', data.plot_field, absolute_frame));

    ax2 = subplot(1,2,2);
    plot_state_field(ax2, X_inside_pred(:,k), data, raw_clim, ...
        sprintf('Full WSINDy recreation %s, Frame %d', data.plot_field, absolute_frame));

    drawnow;
end

end

function clim_vals = robust_raw_clim(X_raw_true, data)

all_values = [];

for k = 1:size(X_raw_true, 2)
    vals = state_to_field(X_raw_true(:,k), data);
    vals = real(vals(:));
    vals = vals(isfinite(vals));
    all_values = [all_values; vals]; %#ok<AGROW>
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
