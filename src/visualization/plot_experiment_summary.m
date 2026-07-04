function plot_experiment_summary(details, data, cfg, k)
%PLOT_EXPERIMENT_SUMMARY Make the standard post-run diagnostic plots.

if nargin < 4 || isempty(k)
    k = 1;
end

plot_l3_error(details, cfg);
plot_full_error(details, cfg);
plot_l3_frame(details, data, k);
plot_full_frame(details, data, k);

end

