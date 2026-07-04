function animate_experiment_summary(details, data, cfg, mrdmd, mode_groups)
%ANIMATE_EXPERIMENT_SUMMARY Run the standard animation set.

if nargin < 5
    mode_groups = default_mode_groups(data);
end

animate_l3_wsindy_comparison(details, data, cfg);
animate_full_wsindy_comparison(details, data, cfg);
animate_mrdmd_mode_groups(mrdmd, data, mode_groups);

end
