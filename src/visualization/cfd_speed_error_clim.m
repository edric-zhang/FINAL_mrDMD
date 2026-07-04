function clim_vals = cfd_speed_error_clim(A, B, npoints)
%CFD_SPEED_ERROR_CLIM Finds reasonable [0, max_error] color scale.

max_val = 0;
for k = 1:min(size(A,2), size(B,2))
    max_val = max(max_val, max(cfd_speed_error(A(:,k), B(:,k), npoints)));
end
if ~isfinite(max_val) || max_val <= 0
    max_val = 1;
end
clim_vals = [0 max_val];

end

