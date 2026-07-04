function speed = cfd_speed_from_state(state, npoints)
%CFD_SPEED_FROM_STATE Convert stacked [u; v] state to speed.

state = real(state(:));
u = state(1:npoints);
v = state(npoints+1:2*npoints);
speed = sqrt(u.^2 + v.^2);

end

