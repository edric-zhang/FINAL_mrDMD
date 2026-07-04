function coef = prune_equation_terms(coef, labels, target_label, relative_tol, max_terms, max_nonlinear_terms)
%PRUNE_EQUATION_TERMS Remove weak and excessive WSINDy equation terms.

if isempty(coef) || all(coef == 0)
    return;
end

abs_coef = abs(coef);
max_coef = max(abs_coef);

if max_coef == 0
    return;
end

small_terms = abs_coef < relative_tol * max_coef;
coef(small_terms) = 0;

nonlinear = false(size(coef));
for ii = 1:length(labels)
    nonlinear(ii) = contains(labels{ii}, '^2') || contains(labels{ii}, ')*(');
end

active_nonlinear = find(coef ~= 0 & nonlinear);
if length(active_nonlinear) > max_nonlinear_terms
    [~, ord] = sort(abs(coef(active_nonlinear)), 'descend');
    drop_idx = active_nonlinear(ord(max_nonlinear_terms+1:end));
    coef(drop_idx) = 0;
end

active = find(coef ~= 0);
if length(active) > max_terms
    [~, ord] = sort(abs(coef(active)), 'descend');
    keep = active(ord(1:max_terms));

    target_self = find(strcmp(labels, target_label), 1);
    if ~isempty(target_self) && coef(target_self) ~= 0 && ~ismember(target_self, keep)
        weakest_keep = keep(end);
        keep(end) = target_self;
        fprintf('Keeping target self term for %s; dropping weaker term %.4e.\n', ...
            target_label, coef(weakest_keep));
    end

    drop = setdiff(active, keep);
    coef(drop) = 0;
end

end

