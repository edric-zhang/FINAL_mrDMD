function animate_l3_wsindy_comparison(details, data, cfg)
%ANIMATE_L3_WSINDY_COMPARISON Animate true vs WSINDy L3 contribution.

plot_l3_wsindy_comparison(details.X_l3_true, details.X_l3_wsindy, ...
    data, ...
    cfg.frames.test_start_idx, size(details.X_l3_true, 2));

end
