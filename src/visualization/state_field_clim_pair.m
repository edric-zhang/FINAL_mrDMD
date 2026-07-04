function clim_vals = state_field_clim_pair(A, B, data, field_type)
%STATE_FIELD_CLIM_PAIR Shared color limits for true vs predicted fields.

if nargin < 4 || isempty(field_type)
    field_type = data.plot_field;
end

max_val = 0;
min_val = inf;

for k = 1:size(A,2)
    values = state_to_field(A(:,k), data, field_type);
    max_val = max(max_val, max(values));
    min_val = min(min_val, min(values));
end
for k = 1:size(B,2)
    values = state_to_field(B(:,k), data, field_type);
    max_val = max(max_val, max(values));
    min_val = min(min_val, min(values));
end

if ~isfinite(max_val) || ~isfinite(min_val) || max_val == min_val
    if strcmpi(field_type, 'speed')
        clim_vals = [0 1];
    else
        clim_vals = [-1 1];
    end
elseif strcmpi(field_type, 'speed') || (strcmpi(field_type, 'scalar') && min_val >= 0)
    clim_vals = [0 max_val];
else
    bound = max(abs([min_val max_val]));
    clim_vals = [-bound bound];
end

end
