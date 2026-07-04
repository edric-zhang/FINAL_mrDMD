function animate_mrdmd_mode_groups(mrdmd, data, mode_groups)
%ANIMATE_MRDMD_MODE_GROUPS Animate selected grouped mrDMD mode contributions.

if nargin < 3 || isempty(mode_groups)
    mode_groups = default_mode_groups(data);
end

plot_cfd_mrdmd_group_animation(mrdmd.list_modes, mrdmd.list_w, mrdmd.list_b, ...
    mrdmd.list_t_start, mrdmd.list_bin_widths, data, mode_groups);

end
