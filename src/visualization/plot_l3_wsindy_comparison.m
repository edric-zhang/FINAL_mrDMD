function plot_l3_wsindy_comparison(X_l3_true, X_l3_wsindy, data, ...
    test_start_idx, recreate_steps)
%PLOT_L3_WSINDY_COMPARISON Animate true vs WSINDy L3 contribution.

figure('Name', 'True vs WSINDy L3 Contribution', ...
    'Position', [100 100 1500 500]);

l3_clim = state_field_clim_pair(X_l3_true, X_l3_wsindy, data);
l3_err_clim = state_field_error_clim(X_l3_true, X_l3_wsindy, data);

for k = 1:recreate_steps
    absolute_frame = test_start_idx + k - 1;

    ax1 = subplot(1,3,1);
    plot_state_field(ax1, X_l3_true(:,k), data, l3_clim, ...
        sprintf('True L3 %s, Frame %d', data.plot_field, absolute_frame));

    ax2 = subplot(1,3,2);
    plot_state_field(ax2, X_l3_wsindy(:,k), data, l3_clim, ...
        sprintf('One-Step WSINDy L3 %s, Frame %d', data.plot_field, absolute_frame));

    ax3 = subplot(1,3,3);
    l3_err_frame = state_field_error(X_l3_true(:,k), X_l3_wsindy(:,k), data);
    cfd_plot_scalar_field(ax3, data.x, data.y, l3_err_frame, l3_err_clim, ...
        sprintf('L3 %s Error, Frame %d', data.plot_field, absolute_frame));

    drawnow;
end

end
