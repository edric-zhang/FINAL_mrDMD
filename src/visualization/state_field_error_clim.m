function clim_vals = state_field_error_clim(A, B, data, field_type)
%STATE_FIELD_ERROR_CLIM Reasonable [0, max_error] color scale.

if nargin < 4 || isempty(field_type)
    field_type = data.plot_field;
end

max_val = 0;
for k = 1:min(size(A,2), size(B,2))
    err = state_field_error(A(:,k), B(:,k), data, field_type);
    max_val = max(max_val, max(err));
end

if ~isfinite(max_val) || max_val <= 0
    max_val = 1;
end
clim_vals = [0 max_val];

end

