function plot_l3_error(details, cfg)
%PLOT_L3_ERROR Plot inside-bin WSINDy recreation error for Level 3 only.

frames = cfg.frames.test_start_idx:cfg.frames.test_end_idx;

figure;
plot(frames, details.l3_recreate_error * 100, 'LineWidth', 1.5);
grid on;
xlabel('Snapshot');
ylabel('L3 Recreation Error %');
title('Inside-Bin WSINDy Recreation Error for Level 3 Only');

end

