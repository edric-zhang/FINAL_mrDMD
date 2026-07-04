function mrdmd = compute_mrdmd(X, dt, cfg)
%COMPUTE_MRDMD Compute multiresolution DMD on snapshot matrix X.

[n, m] = size(X);
L = cfg.mrdmd.L;
maxJ = 2^(L-1);
matrices = cell(1, maxJ);       % List of matrices that we update to run DMD through
matrices{1} = X;




% Tracking variables
list_ml = zeros(L, maxJ);           % List of truncation number for each level, bin
list_b = cell(L, maxJ);             % List of starting vector coefficients for each level, bin mode
list_w = cell(L, maxJ);             % List of eigenvalues for each level, bin mode
list_modes = cell(L, maxJ);         % List of modes for each level, bin
list_bin_widths = zeros(L, maxJ);   % List of bin widths for each level, bin
list_t_start = zeros(L, maxJ);      % List of bin start points for each level, bin
list_t_start(1,1) = 1;
list_bin_widths(1,1) = m;
list_anchor_idx = ones(L, maxJ);      % First meaningful snapshot used to initialize modal amplitudes
if isfield(cfg.mrdmd, 'freq_threshold_cycles_per_snapshot')
    freq_threshold = cfg.mrdmd.freq_threshold_cycles_per_snapshot;
elseif isfield(cfg.mrdmd, 'freq_threshold_hz')
    freq_threshold = cfg.mrdmd.freq_threshold_hz * dt;
else
    freq_threshold = 0.5;
end

for i = 1:L
    J = 2^(i-1);
    next_matrices = cell(1, 2*J);   % Setting an empty list of matrices for the NEXT level

    for j = 1:J
        A = matrices{j};            % Get the next matrix from the previous matrices list
        

        if isempty(A) || size(A, 2) < cfg.mrdmd.min_snapshots_per_bin || any(isnan(A(:)))    % If A is empty or too small,
            next_matrices{2*j - 1} = [];                        % then add empty matrices for the
            next_matrices{2*j} = [];                            % next list and skip the process
            continue;
        end

        current_start = list_t_start(i, j);                     % For (1,1), starts at 0
        current_width = size(A, 2);
        list_bin_widths(i, j) = current_width;                  % For (1,1), starts at 0
        % FINDING THE FIRST NON ZERO FRAME (ANCHOR)
        snapshot_norms = vecnorm(A);
        anchor_tol = 1e-10 * max(snapshot_norms);
        anchor_idx = find(snapshot_norms > anchor_tol, 1, 'first');
        if isempty(anchor_idx)
            next_matrices{2*j - 1} = zeros(n, floor(current_width / 2), 'like', X);
            next_matrices{2*j} = zeros(n, current_width - floor(current_width / 2), 'like', X);
            list_anchor_idx(i, j) = 1;
            continue;
        end
        list_anchor_idx(i, j) = anchor_idx;

        [modes, D, b] = dmd(A, cfg.mrdmd.svd_rank);             % Run DMD
        timeperiod = size(D, 1);                                % Size of D (depends on our r value)

        % Calculate discrete frequencies in cycles per snapshot.
        freqs = zeros(timeperiod, 1);
        for modenum = 1:timeperiod
            lambda = D(modenum, modenum);                       % For each number in our number of modes
            freqs(modenum) = abs(angle(lambda)) / (2 * pi);     % Frequency relative to the saved snapshots
        end

        [sorted_freqs, sort_idx] = sort(freqs);                 % Sort frequencies ascending
        D = D(sort_idx, sort_idx);                              % Sort D, modes, and b in the same way
        modes = modes(:, sort_idx);
        b = single(pinv(modes) * A(:, anchor_idx));

        ml = find(sorted_freqs >= freq_threshold, 1, 'first');  % Find the first frequency to be higher than the threshold

        if isempty(ml)                                          % If none are higher, it flags ALL MODES
            ml = timeperiod + 1;                                % as being SLOW MODES
            slow_inds = 1:timeperiod;
        else
            slow_inds = 1:(ml-1);                               % Or else it flags all up to that index
        end

        if ~isempty(slow_inds)                                  % As long as there exist slow modes
            eigs_slow = diag(D(slow_inds, slow_inds));
            mag = abs(eigs_slow);
            over = mag > 1.0;                                   % If the magnitude of the discrete eigenvalues
            eigs_slow(over) = eigs_slow(over) ./ mag(over);     % is greater than 1, project it onto 1 unit circle

            slowmatrix = zeros(n, current_width, 'like', X);

            time_powers = eigs_slow .^ (0:(current_width - anchor_idx));  % Basically the omega matrix without b
            slowmatrix(:, anchor_idx:end) = modes(:, slow_inds) * ...     % Modes x Omega Matrix
                (b(slow_inds) .* time_powers);
        else
            slowmatrix = zeros(n, current_width, 'like', X);    % If there are no slow modes, just make zeros
        end

        fastmatrix = A - slowmatrix;

        % Storing all the variables
        list_ml(i, j) = ml;
        list_w{i, j} = diag(D(slow_inds, slow_inds));           % Store eigenvalues as a vector
        list_b{i, j} = b(slow_inds);
        list_modes{i, j} = modes(:, slow_inds);

        % Splitting Step: Split the new fastmatrix into half
        midpoint = floor(current_width / 2);
        A1 = fastmatrix(:, 1:midpoint);
        A2 = fastmatrix(:, (midpoint + 1):end);
        left_child_idx = 2*j - 1;
        next_matrices{left_child_idx} = A1;
        right_child_idx = 2*j;
        next_matrices{right_child_idx} = A2;
        if i < L
            list_t_start(i+1, left_child_idx) = current_start;  % Setting starts for the NEXT LEVEL
            list_bin_widths(i+1, left_child_idx) = midpoint;

            list_t_start(i+1, right_child_idx) = current_start + midpoint;
            list_bin_widths(i+1, right_child_idx) = current_width - midpoint;
        end
    end
    matrices = next_matrices; % Update the new matrix list
end

mrdmd.list_ml = list_ml;
mrdmd.list_b = list_b;
mrdmd.list_w = list_w;
mrdmd.list_modes = list_modes;
mrdmd.list_bin_widths = list_bin_widths;
mrdmd.list_t_start = list_t_start;
mrdmd.L = L;
mrdmd.maxJ = maxJ;
mrdmd.freq_threshold_cycles_per_snapshot = freq_threshold;
mrdmd.list_anchor_idx = list_anchor_idx;

end
