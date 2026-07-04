function plot_cfd_mrdmd_group_animation(list_modes, list_w, list_b, list_t_start, list_bin_widths, ...
    data, mode_groups)
%PLOT_CFD_MRDMD_GROUP_ANIMATION Animate chosen grouped mrDMD mode contributions.

num_groups = size(mode_groups, 1);
num_cols = 2;
num_rows = ceil(num_groups / num_cols);

group_frames = cell(num_groups, 1);
valid_groups = false(num_groups, 1);
max_steps = 0;

for gg = 1:num_groups
    specs = mode_groups{gg, 1};
    frames = [];

    for ii = 1:size(specs, 1)
        lev = specs(ii, 1);
        bin = specs(ii, 2);
        mode_idx = specs(ii, 3);

        if ~is_valid_mrdmd_mode(list_modes, list_w, lev, bin, mode_idx)
            continue;
        end

        t_start = list_t_start(lev, bin);
        bin_width = list_bin_widths(lev, bin);
        t_end = t_start + bin_width - 1;
        frames = union(frames, t_start:t_end);
    end

    if ~isempty(frames)
        group_frames{gg} = frames;
        valid_groups(gg) = true;
        max_steps = max(max_steps, length(frames));
    end
end

clim_by_group = cfd_group_animation_clim(list_modes, list_w, list_b, ...
    list_t_start, list_bin_widths, data, mode_groups, group_frames, valid_groups);

figure('Name', 'Grouped mrDMD Mode Contributions', ...
    'Position', [100 30 1500 450*num_rows]);

tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

scatter_handles = gobjects(num_groups, 1);
title_handles = gobjects(num_groups, 1);
axis_handles = gobjects(num_groups, 1);

for gg = 1:num_groups
    ax = nexttile;
    axis_handles(gg) = ax;

    field_type = lower(mode_groups{gg, 2});
    group_title = mode_groups{gg, 3};

    values = nan(data.npoints, 1);
    if valid_groups(gg)
        state = cfd_group_contribution_state(list_modes, list_w, list_b, ...
            list_t_start, list_bin_widths, data.n, mode_groups{gg, 1}, group_frames{gg}(1));
        values = state_to_field(state, data, field_type);
    end

    scatter_handles(gg) = scatter(ax, data.x, data.y, 12, values, 'filled');
    axis(ax, 'equal');
    axis(ax, 'tight');
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    colorbar(ax);
    set(ax, 'CLim', clim_by_group(gg, :));
    title_handles(gg) = title(ax, sprintf('[%s] %s', field_type, group_title));
end

for step = 1:max_steps
    for gg = 1:num_groups
        field_type = lower(mode_groups{gg, 2});
        group_title = mode_groups{gg, 3};

        if ~valid_groups(gg)
            set(scatter_handles(gg), 'CData', nan(data.npoints, 1));
            set(title_handles(gg), 'String', sprintf('[%s] %s | missing', field_type, group_title));
            continue;
        end

        frames = group_frames{gg};
        frame = frames(min(step, length(frames)));

        state = cfd_group_contribution_state(list_modes, list_w, list_b, ...
            list_t_start, list_bin_widths, data.n, mode_groups{gg, 1}, frame);
        values = state_to_field(state, data, field_type);

        set(scatter_handles(gg), 'CData', values);
        set(axis_handles(gg), 'CLim', clim_by_group(gg, :));
        set(title_handles(gg), 'String', sprintf('[%s] %s | Frame %d/%d-%d', ...
            field_type, group_title, frame, frames(1), frames(end)));
    end

    drawnow;
end

end
