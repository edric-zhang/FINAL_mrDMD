function available_modes = print_available_modes(mrdmd, plot_start_idx, plot_end_idx)
%PRINT_AVAILABLE_MODES Print mode coordinates overlapping a frame range.

fprintf('\n===========================================================');
fprintf('\n   AVAILABLE COORDINATES FOR FRAMES %d TO %d', plot_start_idx, plot_end_idx);
fprintf('\n===========================================================\n');

available_modes = cell(mrdmd.L, mrdmd.maxJ);

for i = 1:mrdmd.L
    J = 2^(i-1);
    for j = 1:J
        if isempty(mrdmd.list_modes{i, j}), continue; end

        t_start = mrdmd.list_t_start(i, j);
        bin_width = mrdmd.list_bin_widths(i, j);
        t_end = t_start + bin_width - 1;

        % Check if this bin overlaps with our visual playback window
        if t_start <= plot_end_idx && t_end >= plot_start_idx
            num_modes_available = size(mrdmd.list_modes{i, j}, 2);
            available_modes{i,j} = mrdmd.list_modes{i,j};
            fprintf('Level %d, Bin %d (Frames %d to %d)\n', i, j, t_start, t_end);
            fprintf('  Available Mode Indices: [');
            for m_idx = 1:num_modes_available
                if m_idx == num_modes_available
                    fprintf('%d', m_idx);
                else
                    fprintf('%d, ', m_idx);
                end
            end
            fprintf(']\n\n');
        end
    end
end
fprintf('===========================================================\n\n');

end

