function clim_vals = cfd_speed_clim_pair(A, B, npoints)
%CFD_SPEED_CLIM_PAIR Shared color limit identifier for true vs predicted.

max_val = 0;
for k = 1:size(A,2)
    max_val = max(max_val, max(cfd_speed_from_state(A(:,k), npoints)));
end
for k = 1:size(B,2)
    max_val = max(max_val, max(cfd_speed_from_state(B(:,k), npoints)));
end
if ~isfinite(max_val) || max_val <= 0
    max_val = 1;
end
clim_vals = [0 max_val];

end

