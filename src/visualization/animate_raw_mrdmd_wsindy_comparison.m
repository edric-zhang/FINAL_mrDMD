function animate_raw_mrdmd_wsindy_comparison(details, data, cfg)
%ANIMATE_RAW_MRDMD_WSINDY_COMPARISON Animate raw, true mrDMD, and WSINDy recreation.

if ~isfield(details, 'X_raw_true')
    error('details.X_raw_true is missing. Rerun the experiment to regenerate details.');
end

plot_raw_mrdmd_wsindy_comparison(details.X_raw_true, details.X_inside_true, ...
    details.X_inside_pred, data, cfg.frames.test_start_idx, size(details.X_raw_true, 2));

end
