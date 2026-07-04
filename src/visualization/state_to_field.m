function values = state_to_field(state, data, field_type)
%STATE_TO_FIELD Convert one state vector to plottable scalar values.
%
% For uv data:
%   field_type = 'u', 'v', or 'speed'
% For scalar data:
%   field_type = 'scalar'

if nargin < 3 || isempty(field_type)
    field_type = data.plot_field;
end

field_type = lower(field_type);
layout = lower(data.state_layout);
state = real(state(:));

switch layout
    case 'uv'
        u = state(1:data.npoints);
        v = state(data.npoints+1:2*data.npoints);

        switch field_type
            case 'u'
                values = u;
            case 'v'
                values = v;
            case 'speed'
                values = sqrt(u.^2 + v.^2);
            otherwise
                error('For uv data, field_type must be ''u'', ''v'', or ''speed''.');
        end

    case 'scalar'
        if ~strcmp(field_type, 'scalar')
            error('For scalar data, field_type must be ''scalar''.');
        end
        values = state(1:data.npoints);

    otherwise
        error('Unknown data.state_layout "%s".', data.state_layout);
end

end

