function plot_cancellation_metrics(details, cfg)
%PLOT_CANCELLATION_METRICS Show whether the target is a small part or a cancellation partner.

frames = cfg.frames.test_start_idx + (0:size(details.X_raw_true, 2)-1);
background = details.X_inside_true - details.X_l3_true;

raw_norm = vecnorm(details.X_raw_true);
mrdmd_norm = vecnorm(details.X_inside_true);
target_norm = vecnorm(details.X_l3_true);
background_norm = vecnorm(background);
true_full_norm = vecnorm(background + details.X_l3_true);
pred_full_norm = vecnorm(background + details.X_l3_wsindy);

figure('Name', 'Cancellation Metrics', 'NumberTitle', 'off');
plot(frames, raw_norm, 'LineWidth', 1.5, 'DisplayName', 'raw');
hold on;
plot(frames, mrdmd_norm, '--', 'LineWidth', 1.2, 'DisplayName', 'mrDMD full');
plot(frames, target_norm, 'LineWidth', 1.4, 'DisplayName', sprintf('target L%d', cfg.wsindy.target_level));
plot(frames, background_norm, 'LineWidth', 1.4, 'DisplayName', 'background');
plot(frames, true_full_norm, ':', 'LineWidth', 1.5, 'DisplayName', 'background + true target');
plot(frames, pred_full_norm, 'LineWidth', 1.8, 'DisplayName', 'background + WSINDy target');
grid on;
xlabel('Snapshot');
ylabel('Frame norm');
title('Cancellation Balance');
legend('Location', 'best');

end
