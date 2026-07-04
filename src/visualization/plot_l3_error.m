function plot_l3_error(details, cfg)
%PLOT_L3_ERROR Plot inside-bin WSINDy recreation error for Level 3 only.

err = details.l3_recreate_error(:).';
frames = cfg.frames.test_start_idx + (0:numel(err)-1);

figure('Name', 'L3 Error', 'NumberTitle', 'off');
plot(frames, err * 100, 'LineWidth', 1.5, 'DisplayName', 'WSINDy vs selected L3');
hold on;

if isfield(details, 'raw_l3_error')
    raw_err = details.raw_l3_error(:).';
    raw_frames = cfg.frames.test_start_idx + (0:numel(raw_err)-1);
    plot(raw_frames, raw_err * 100, 'LineWidth', 1.5, ...
        'DisplayName', 'L3-only WSINDy vs raw (diagnostic)');
end

grid on;
xlabel('Snapshot');
ylabel('L3 Error (% of Median Snapshot Norm)');
title('Inside-Bin WSINDy Recreation Error for Level 3 Only');
legend('Location', 'best');

end
