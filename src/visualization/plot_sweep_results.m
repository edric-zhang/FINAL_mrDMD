function plot_sweep_results(result_table)
%PLOT_SWEEP_RESULTS Plot standard WSINDy sweep diagnostics.

figure;
tiledlayout(1,2);

ax1 = nexttile;
plot(ax1, result_table.mean_L3_error, 'o-', 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 'Sweep Run');
ylabel(ax1, 'Mean L3 Error %');
title(ax1, 'Sweep Mean L3 Error');

ax2 = nexttile;
plot(ax2, result_table.mean_full_error, 'o-', 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel(ax2, 'Sweep Run');
ylabel(ax2, 'Mean Full Error %');
title(ax2, 'Sweep Mean Full CFD Error');

end

