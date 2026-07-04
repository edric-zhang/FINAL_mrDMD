function plot_full_error(details, cfg)
%PLOT_FULL_ERROR Plot inside-bin full CFD recreation error.

frames = cfg.frames.test_start_idx:cfg.frames.test_end_idx;

figure;
plot(frames, details.inside_full_error * 100, 'LineWidth', 1.5);
grid on;
xlabel('Snapshot');
ylabel('Full CFD Recreation Error %');
title('Inside-Bin Full CFD Recreation Error');

end

