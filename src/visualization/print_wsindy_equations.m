function print_wsindy_equations(w_second, labels_second, mode_labels, target_cols)
%PRINT_WSINDY_EQUATIONS Print discovered Level 3 derivative equations.

for kk = 1:length(target_cols)

    target_col = target_cols(kk);
    coef_vector = w_second{kk};
    labels2 = labels_second{kk};

    print_tol = 1e-6 * max(abs(coef_vector));
    active_idx = find(abs(coef_vector) > print_tol);

    fprintf('\nEquation for d/dt of %s:\n', mode_labels{target_col});

    if isempty(active_idx)
        fprintf('    d(%s)/dt = 0\n', mode_labels{target_col});
        continue;
    end

    equation_str = sprintf('    d(%s)/dt = ', mode_labels{target_col});

    for q = 1:length(active_idx)
        row = active_idx(q);
        weight = coef_vector(row);

        if row <= length(labels2)
            label = labels2{row};
        else
            label = sprintf('Theta_%d', row);
        end

        if q == 1
            equation_str = sprintf('%s %.4e * %s', equation_str, weight, label);
        else
            if weight >= 0
                equation_str = sprintf('%s + %.4e * %s', equation_str, weight, label);
            else
                equation_str = sprintf('%s - %.4e * %s', equation_str, abs(weight), label);
            end
        end
    end

    fprintf('%s\n', equation_str);
end

fprintf('\n==========================================================================\n');

end

