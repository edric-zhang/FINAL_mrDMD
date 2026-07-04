function state = cfd_group_contribution_state(list_modes, list_w, list_b, ...
    list_t_start, list_bin_widths, state_size, mode_specs, frame)
%CFD_GROUP_CONTRIBUTION_STATE Build the time evolution of chosen mrDMD modes.
% Use Mode * b * lambda^time.

state = zeros(state_size, 1);

for ii = 1:size(mode_specs, 1)
    lev = mode_specs(ii, 1);
    bin = mode_specs(ii, 2);
    mode_idx = mode_specs(ii, 3);

    if ~is_valid_mrdmd_mode(list_modes, list_w, lev, bin, mode_idx)
        continue;
    end

    t_start = list_t_start(lev, bin);
    bin_width = list_bin_widths(lev, bin);
    t_end = t_start + bin_width - 1;

    if frame < t_start || frame > t_end
        continue;
    end

    rel_time = frame - t_start;
    amp = list_b{lev, bin}(mode_idx) * list_w{lev, bin}(mode_idx)^rel_time;
    state = state + real(list_modes{lev, bin}(:, mode_idx) * amp);
end

end
