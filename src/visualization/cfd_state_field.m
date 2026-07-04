function values = cfd_state_field(state, npoints, field_type)
%CFD_STATE_FIELD Convert stacked state to 'u', 'v', or 'speed' values.

u = real(state(1:npoints));
v = real(state(npoints+1:2*npoints));

switch lower(field_type)
    case 'u'
        values = u;
    case 'v'
        values = v;
    case 'speed'
        values = sqrt(u.^2 + v.^2);
    otherwise
        error('field_type must be ''u'', ''v'', or ''speed''.');
end

end

