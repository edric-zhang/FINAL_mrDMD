function err = state_field_error(state_true, state_pred, data, field_type)
%STATE_FIELD_ERROR Pointwise error in the selected plottable field.

if nargin < 4 || isempty(field_type)
    field_type = data.plot_field;
end

true_values = state_to_field(state_true, data, field_type);
pred_values = state_to_field(state_pred, data, field_type);
err = abs(true_values - pred_values);

end
