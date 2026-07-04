function cfd_plot_state_speed(ax, state, x, y, npoints, clim_vals, title_text)
%CFD_PLOT_STATE_SPEED Converts [u;v] state to speed, then plots.

speed = cfd_speed_from_state(state, npoints);
cfd_plot_scalar_field(ax, x, y, speed, clim_vals, title_text);

end

