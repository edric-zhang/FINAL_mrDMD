function dydt = mrdmd_wsindy_inside_rhs(tt, y, w_second, labels_second, ...
    mode_labels, target_cols, xobs, tobs)
%MRDMD_WSINDY_INSIDE_RHS Propagate target L3 modal amplitudes.

dydt = zeros(length(target_cols),1);

for kk = 1:length(target_cols)

    coef = w_second{kk};
    labels = labels_second{kk};

    theta = zeros(length(labels),1);

    for r = 1:length(labels)
        theta(r) = evaluate_inside_label(labels{r}, tt, y, ...
            mode_labels, target_cols, xobs, tobs);
    end

    dydt(kk) = coef(:).' * theta(:);
end

end

