function plot_full_metrics(details, cfg)
%PLOT_FULL_METRICS Plot whole-frame WSINDy quality against raw data.

frames = cfg.frames.test_start_idx + (0:size(details.X_raw_true, 2)-1);

mrdmd_frame_error = zeros(size(frames));
for k = 1:numel(frames)
    mrdmd_frame_error(k) = norm(details.X_raw_true(:,k) - details.X_inside_true(:,k)) / ...
        (norm(details.X_raw_true(:,k)) + eps);
end

figure('Name', 'Full Metrics', 'NumberTitle', 'off');

subplot(2,1,1);
plot(frames, details.raw_full_frame_error * 100, 'LineWidth', 1.6, ...
    'DisplayName', 'WSINDy vs raw');
hold on;
plot(frames, mrdmd_frame_error * 100, '--', 'LineWidth', 1.3, ...
    'DisplayName', 'mrDMD vs raw');
grid on;
xlabel('Snapshot');
ylabel('Error (% raw norm)');
title('Full-Orchestra Error');
legend('Location', 'best');

subplot(2,1,2);
plot(frames, details.raw_full_cos, 'LineWidth', 1.6);
grid on;
xlabel('Snapshot');
ylabel('Cosine');
ylim([-1 1]);
title('Full-Orchestra Direction Agreement');

end
