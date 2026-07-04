function plot_full_wsindy_comparison(X_inside_true, X_inside_pred, data, ...
    test_start_idx, recreate_steps)
%PLOT_FULL_WSINDY_COMPARISON Animate true vs predicted full inside-bin state.

figure('WindowStyle', 'normal', ...
    'Name', 'Full Inside-Bin Reconstruction', ...
    'NumberTitle', 'off', ...
    'Position', [100 100 1200 500]);

% Use the true/original field to define the visual scale so prediction
% outliers do not wash out the original panel.
full_clim = robust_true_clim(X_inside_true, data);

for k = 1:recreate_steps
    absolute_frame = test_start_idx + k - 1;

    ax1 = subplot(1,2,1);
    plot_state_field(ax1, X_inside_true(:,k), data, full_clim, ...
        sprintf('Full mrDMD %s, Frame %d', data.plot_field, absolute_frame));

    ax2 = subplot(1,2,2);
    plot_state_field(ax2, X_inside_pred(:,k), data, full_clim, ...
        sprintf('Inside-Bin WSINDy %s, Frame %d', data.plot_field, absolute_frame));

    drawnow;
end

end

function clim_vals = robust_true_clim(X_true, data)

all_values = [];

for k = 1:size(X_true, 2)
    vals = state_to_field(X_true(:,k), data);
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
