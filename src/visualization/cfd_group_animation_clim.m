function clim_by_group = cfd_group_animation_clim(list_modes, list_w, list_b, ...
    list_t_start, list_bin_widths, data, mode_groups, group_frames, valid_groups)
%CFD_GROUP_ANIMATION_CLIM Compute fixed color limits for mode animation panels.

num_groups = size(mode_groups, 1);
clim_by_group = repmat([-1 1], num_groups, 1);

for gg = 1:num_groups
    if ~valid_groups(gg)
        continue;
    end

    field_type = lower(mode_groups{gg, 2});
    max_abs_val = 0;
    frames = group_frames{gg};

    for kk = 1:length(frames)
        state = cfd_group_contribution_state(list_modes, list_w, list_b, ...
            list_t_start, list_bin_widths, data.n, mode_groups{gg, 1}, frames(kk));
        values = state_to_field(state, data, field_type);
        max_abs_val = max(max_abs_val, max(abs(values)));
    end

    if ~isfinite(max_abs_val) || max_abs_val <= 0
        max_abs_val = 1;
    end

    if strcmpi(field_type, 'speed')
        clim_by_group(gg, :) = [0 max_abs_val];
    else
        clim_by_group(gg, :) = [-max_abs_val max_abs_val];
    end
end

end
