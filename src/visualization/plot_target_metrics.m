function plot_target_metrics(details, cfg)
%PLOT_TARGET_METRICS Plot selected target-level replacement diagnostics.

frames = cfg.frames.test_start_idx + (0:size(details.X_l3_true, 2)-1);
target_name = sprintf('Target L%d', cfg.wsindy.target_level);
target_norm = vecnorm(details.X_l3_true);
pred_norm = vecnorm(details.X_l3_wsindy);

figure('Name', 'Target Metrics', 'NumberTitle', 'off');

subplot(3,1,1);
plot(frames, details.l3_recreate_error * 100, 'LineWidth', 1.5);
grid on;
xlabel('Snapshot');
ylabel('Error (% global)');
title([target_name ' Error vs mrDMD']);

subplot(3,1,2);
plot(frames, details.l3_spatial_corr, 'LineWidth', 1.5);
grid on;
xlabel('Snapshot');
ylabel('Spatial corr');
ylim([-1 1]);
title([target_name ' Direction Agreement']);

subplot(3,1,3);
plot(frames, target_norm, 'LineWidth', 1.5, 'DisplayName', 'mrDMD target');
hold on;
plot(frames, pred_norm, '--', 'LineWidth', 1.5, 'DisplayName', 'WSINDy target');
grid on;
xlabel('Snapshot');
ylabel('Norm');
title([target_name ' Loudness']);
legend('Location', 'best');

end
