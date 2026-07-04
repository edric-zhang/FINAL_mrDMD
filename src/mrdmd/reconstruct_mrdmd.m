function [X_rec, mrdmd_error] = reconstruct_mrdmd(mrdmd, X, frame_start, numsnapshots)
%Reconstruct the mrDMD contribution and snapshot errors.
    
    [n, m] = size(X);
    X_rec = zeros(n, m, 'like', X);
    
    for i = 1:mrdmd.L
        J = 2^(i-1);
        for j = 1:J
            if isempty(mrdmd.list_modes{i, j})
                continue;
            end
            modes = mrdmd.list_modes{i, j};
            eigs_slow = mrdmd.list_w{i, j};
            b = mrdmd.list_b{i, j};
    
            t_start = mrdmd.list_t_start(i, j);                     % Setting start/end - modes only exist locally
            bin_width = mrdmd.list_bin_widths(i, j);
    
            if bin_width == 0
                continue;
            end
    
            t_end = t_start + bin_width - 1;
            mag = abs(eigs_slow);
            over = mag > 1.0;                                       % Project eigs>1 to 1
            eigs_slow(over) = eigs_slow(over) ./ mag(over);
    
            anchor_idx = 1;
            if isfield(mrdmd, 'list_anchor_idx')
                anchor_idx = mrdmd.list_anchor_idx(i, j);
            end
            
            local_rec = zeros(n, bin_width, 'like', X);
            time_powers = eigs_slow .^ (0:(bin_width - anchor_idx));
            local_rec(:, anchor_idx:end) = modes * (b .* time_powers);
    
            % Protect against matrix index rounding overflows at tree boundaries
            if t_end > m                                            % Make overflow over m into m
                t_end = m;
                local_rec = local_rec(:, 1:(t_end - t_start + 1));
            end
    
            X_rec(:, t_start:t_end) = X_rec(:, t_start:t_end) + local_rec; % Add on the local mode data
        end
    end
    X_rec = real(X_rec);
    
    % Error Graph
    % Error Graph
    mrdmd_error = NaN(1, numsnapshots);
    
    snapshot_norms = vecnorm(X);
    valid_norms = snapshot_norms(snapshot_norms > 0);
    
    % Fixed scale: robust typical earthquake-frame magnitude
    global_scale = median(valid_norms);
    
    for k = 1:m
        true_snapshot = X(:, k);
        rec_snapshot = X_rec(:, k);
    
        absolute_idx = frame_start + k - 1;
    
        mrdmd_error(absolute_idx) = norm(true_snapshot - rec_snapshot) / (global_scale + eps);
    end

end

