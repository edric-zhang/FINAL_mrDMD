function plot_reconstruction_error(mrdmd_error, numsnapshots)
%PLOT_RECONSTRUCTION_ERROR Plot snapshot-wise mrDMD reconstruction error.

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);

% Plotting against the absolute indices preserves the alignment layout
plot(mrdmd_error * 100, 'r-', 'LineWidth', 1.5, 'DisplayName', 'MRDMD');
hold on;
xlim([1, numsnapshots]); % Lock the visual frame boundary to the total absolute snapshots
ylim([0 200]);
grid on;
legend('Location', 'best');
xlabel('Snapshot (Absolute Timeline)');
ylabel('Error %');
title('Snapshot-wise Reconstruction Performance');

end

