function animate_full_wsindy_comparison(details, data, cfg)
%ANIMATE_FULL_WSINDY_COMPARISON Animate true vs predicted full inside-bin CFD state.

plot_full_wsindy_comparison(details.X_inside_true, details.X_inside_pred, ...
    data, ...
    cfg.frames.test_start_idx, size(details.X_inside_true, 2));

end
