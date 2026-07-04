function plot_full_wsindy_comparison(X_inside_true, X_inside_pred, data, ...
    test_start_idx, recreate_steps)
%PLOT_FULL_WSINDY_COMPARISON Animate true vs predicted full inside-bin CFD state.

figure('Name', 'Full Inside-Bin CFD Reconstruction', ...
    'Position', [100 100 1200 500]);

full_clim = state_field_clim_pair(X_inside_true, X_inside_pred, data);

for k = 1:recreate_steps
    absolute_frame = test_start_idx + k - 1;

    ax1 = subplot(1,2,1);
    plot_state_field(ax1, X_inside_true(:,k), data, full_clim, ...
        sprintf('Original %s, Frame %d', data.plot_field, absolute_frame));

    ax2 = subplot(1,2,2);
    plot_state_field(ax2, X_inside_pred(:,k), data, full_clim, ...
        sprintf('Inside-Bin WSINDy %s, Frame %d', data.plot_field, absolute_frame));

    drawnow;
end

end
