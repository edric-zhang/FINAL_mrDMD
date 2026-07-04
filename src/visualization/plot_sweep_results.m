function plot_sweep_results(result_table)
%PLOT_SWEEP_RESULTS Plot standard WSINDy sweep diagnostics.

figure;
tiledlayout(1,3);

ax1 = nexttile;
plot(ax1, result_table.mean_L3_error, 'o-', 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 'Sweep Run');
ylabel(ax1, 'Mean L3 Error %');
title(ax1, 'Sweep Mean L3 Error');

ax2 = nexttile;
plot(ax2, result_table.mean_raw_full_error, 'o-', 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel(ax2, 'Sweep Run');
ylabel(ax2, 'Mean Raw Full Error %');
title(ax2, 'Sweep Mean Raw Full Error');

ax3 = nexttile;
plot(ax3, result_table.mean_raw_full_corr, 'o-', 'LineWidth', 1.5);
grid(ax3, 'on');
xlabel(ax3, 'Sweep Run');
ylabel(ax3, 'Mean Raw Full Corr');
title(ax3, 'Sweep Mean Raw Full Correlation');

end
