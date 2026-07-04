function [xobs, tobs] = build_mrdmd_xobs_for_labels( ...
    list_w, list_b, list_t_start, list_bin_widths, ...
    dt, ~, mode_labels, start_idx, end_idx)
%BUILD_MRDMD_XOBS_FOR_LABELS Builds modal amplitudes based on b * lambda^t.

m_interval = end_idx - start_idx + 1;
xobs = zeros(m_interval, length(mode_labels));

for c = 1:length(mode_labels)

    [lev, bin, mode_idx] = parse_mode_label(mode_labels{c});

    eigs_slow = list_w{lev, bin};
    b = list_b{lev, bin};
    bin_width = list_bin_widths(lev, bin);
    t_start = list_t_start(lev, bin);
    t_end = t_start + bin_width - 1;

    if start_idx > t_end || end_idx < t_start
        continue;
    end

    local_start = max(start_idx, t_start);
    local_end = min(end_idx, t_end);

    rel_start = local_start - t_start;
    rel_end = local_end - t_start;

    insert_start = local_start - start_idx + 1;
    insert_end = local_end - start_idx + 1;

    a_local = b(mode_idx) * eigs_slow(mode_idx).^(rel_start:rel_end);
    xobs(insert_start:insert_end, c) = real(a_local(:));
end

tobs = (0:m_interval-1)' * dt;

end
