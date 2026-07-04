function err_speed = cfd_speed_error(state_true, state_pred, npoints)
%CFD_SPEED_ERROR Pointwise velocity error from stacked [u; v] states.

state_true = real(state_true(:));
state_pred = real(state_pred(:));
du = state_true(1:npoints) - state_pred(1:npoints);
dv = state_true(npoints+1:2*npoints) - state_pred(npoints+1:2*npoints);
err_speed = sqrt(du.^2 + dv.^2);

end

