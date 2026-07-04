function plot_state_field(ax, state, data, clim_vals, title_text, field_type)
%PLOT_STATE_FIELD Plot the selected field from a state vector.

if nargin < 6 || isempty(field_type)
    field_type = data.plot_field;
end

values = state_to_field(state, data, field_type);
cfd_plot_scalar_field(ax, data.x, data.y, values, clim_vals, title_text);

end
