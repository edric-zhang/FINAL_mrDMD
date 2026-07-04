function plot_full_error(details, cfg)
%PLOT_FULL_ERROR Plot raw and modal full-recreation errors.

figure('Name', 'Full Error', 'NumberTitle', 'off');
hold on;

if isfield(details, 'raw_full_error')
    raw_err = details.raw_full_error(:).';
    raw_frames = cfg.frames.test_start_idx + (0:numel(raw_err)-1);
    plot(raw_frames, raw_err * 100, 'LineWidth', 1.8, ...
        'DisplayName', 'WSINDy full vs raw');
else
    warning('details.raw_full_error is missing. Rerun the experiment to regenerate details.');
end

if isfield(details, 'inside_full_error')
    modal_err = details.inside_full_error(:).';
    modal_frames = cfg.frames.test_start_idx + (0:numel(modal_err)-1);
    plot(modal_frames, modal_err * 100, '--', 'LineWidth', 1.5, ...
        'DisplayName', 'WSINDy full vs true mrDMD');
else
    warning('details.inside_full_error is missing. Rerun the experiment to regenerate details.');
end

grid on;
xlabel('Snapshot');
ylabel('Full Error (% of Median Snapshot Norm)');
title('Inside-Bin Full Recreation Error');
legend('Location', 'best');

end
