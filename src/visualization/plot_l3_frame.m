function plot_l3_frame(details, data, k)
%PLOT_L3_FRAME Plot true and WSINDy-predicted L3 contribution at test-frame index k.

if nargin < 3 || isempty(k)
    k = 1;
end

clim_vals = state_field_clim_pair(details.X_l3_true, details.X_l3_wsindy, data);
err_clim = state_field_error_clim(details.X_l3_true, details.X_l3_wsindy, data);
err_vals = state_field_error(details.X_l3_true(:,k), details.X_l3_wsindy(:,k), data);
if isfield(data, 'colormap')
    cmap = data.colormap;
else
    cmap = [];
end

figure;
tiledlayout(1,3);

ax1 = nexttile;
plot_state_field(ax1, details.X_l3_true(:,k), data, clim_vals, sprintf('True L3 %s', data.plot_field));

ax2 = nexttile;
plot_state_field(ax2, details.X_l3_wsindy(:,k), data, clim_vals, sprintf('WSINDy L3 %s', data.plot_field));

ax3 = nexttile;
cfd_plot_scalar_field(ax3, data.x, data.y, err_vals, err_clim, sprintf('L3 %s Error', data.plot_field), cmap);

end
