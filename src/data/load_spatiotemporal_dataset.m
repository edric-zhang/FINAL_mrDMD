function data = load_spatiotemporal_dataset(cfg)
%LOAD_SPATIOTEMPORAL_DATASET Load spatial-state snapshots over time.

fprintf('Loading %s...\n', cfg.data.mat_file);
raw = load(cfg.data.mat_file);

data = validate_spatiotemporal_dataset(raw, cfg);
data.X = single(data.X);

X_all = data.X;
if size(X_all, 2) >= cfg.data.max_working_snapshots + 1
    data.X_test = single(X_all(:, cfg.data.max_working_snapshots+1:end));
else
    data.X_test = single([]);
end

% Keep the same 1:400 working window used by the original sweep.
% Increase this later once you want to exploit more snapshots.
data.X = X_all(:, 1:min(cfg.data.max_working_snapshots, size(X_all, 2)));

frame_start = cfg.data.frame_start;
frame_end = size(data.X, 2);
data.X = data.X(:, frame_start:frame_end);
data.X_total_size = size(data.X, 2);
data.n = size(data.X, 1);
data.m = size(data.X, 2);
data.numsnapshots = size(data.X, 2);   % keep this as full working snapshot length
data.total_time = data.dt * (data.m - 1);  % True total duration of analyzed snapshots
data.t = (0:data.dt:data.total_time)';

end
