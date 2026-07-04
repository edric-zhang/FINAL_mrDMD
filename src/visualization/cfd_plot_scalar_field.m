function cfd_plot_scalar_field(ax, x, y, values, clim_vals, title_text)
%CFD_PLOT_SCALAR_FIELD Makes scatter plot over (x,y) using scalar values as color.

scatter(ax, x, y, 12, real(values(:)), 'filled');
axis(ax, 'equal');
axis(ax, 'tight');
xlabel(ax, 'x');
ylabel(ax, 'y');
title(ax, title_text);
colorbar(ax);
set(ax, 'CLim', clim_vals);

end

