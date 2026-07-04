function plot_full_frame(details, data, k)
%PLOT_FULL_FRAME Plot true, predicted, and speed error for full inside-bin CFD state.

if nargin < 3 || isempty(k)
    k = 1;
end

clim_vals = state_field_clim_pair(details.X_inside_true, details.X_inside_pred, data);
err_clim = state_field_error_clim(details.X_inside_true, details.X_inside_pred, data);
err_vals = state_field_error(details.X_inside_true(:,k), details.X_inside_pred(:,k), data);
if isfield(data, 'colormap')
    cmap = data.colormap;
else
    cmap = [];
end

figure;
tiledlayout(1,3);

ax1 = nexttile;
plot_state_field(ax1, details.X_inside_true(:,k), data, clim_vals, sprintf('True %s', data.plot_field));

ax2 = nexttile;
plot_state_field(ax2, details.X_inside_pred(:,k), data, clim_vals, sprintf('Predicted %s', data.plot_field));

ax3 = nexttile;
cfd_plot_scalar_field(ax3, data.x, data.y, err_vals, err_clim, sprintf('%s Error', data.plot_field), cmap);

end
