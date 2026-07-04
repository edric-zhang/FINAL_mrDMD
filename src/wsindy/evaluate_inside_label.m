function val = evaluate_inside_label(label, tt, y, mode_labels, target_cols, xobs, tobs)
%EVALUATE_INSIDE_LABEL Turns labels like '(L2 B2 Mode11)^2' into values.
% Handles self quadratics and cross terms for now.

% self quadratic
if startsWith(label, '(') && endsWith(label, ')^2')
    inner = extractBetween(label, "(", ")^2");
    base = evaluate_inside_label(inner{1}, tt, y, mode_labels, target_cols, xobs, tobs);
    val = base^2;
    return;
end

% cross term: (mode A)*(mode B)
if contains(label, ')*(')
    parts = regexp(label, '^\((.*)\)\*\((.*)\)$', 'tokens');
    a = parts{1}{1};
    b = parts{1}{2};

    va = evaluate_inside_label(a, tt, y, mode_labels, target_cols, xobs, tobs);
    vb = evaluate_inside_label(b, tt, y, mode_labels, target_cols, xobs, tobs);

    val = va * vb;
    return;
end

% if label is a target L3 variable, use integrated y
idx_target = find(strcmp(mode_labels(target_cols), label), 1);

if ~isempty(idx_target)
    val = y(idx_target);
    return;
end

% otherwise, this is an input/ancestor mode.
% Inside-bin recreation uses the observed xobs value at time tt.
idx = find(strcmp(mode_labels, label), 1);

if isempty(idx)
    error('Could not find label: %s', label);
end

val = interp1(tobs, xobs(:,idx), tt, 'linear', 'extrap');

end

