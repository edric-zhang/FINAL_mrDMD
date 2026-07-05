function [w_second, labels_second, mode_labels, target_cols, mu_xobs, sigma_xobs, w_linear] = run_mrdmd_wsindy( ...
    list_w, list_b, list_modes, list_t_start, list_bin_widths, ...
    dt, m, input_levels, target_level, plot_start_idx, plot_end_idx, ...
    top_input_modes_per_level, top_target_modes, lambda1, lambda2, gamma, ...
    max_terms_per_equation, max_quadratic_base_terms, list_anchor_idx, include_constant_term, include_time_term, force_target_peer_terms, max_peer_l3_per_equation)
%RUN_MRDMD_WSINDY Fit two-pass WSINDy models on mrDMD modal amplitudes.
%
% Builds modal amplitude timeseries (xobs) from stored DMD modes.
% Keeps only modes in selected interval.
% Ranks and only keeps mode based on RMS Energy.
% Keeps TOP INPUT MODES from L1/L2 and TOP TARGET MODES from L3.
% First pass WSINDY - gets relevant modes - computes more detailed library
% with those modes - runs Second pass WSINDY.
% Caps large coefficients, caps self-growth, prunes weak/excess terms.

if nargin < 14 || isempty(lambda1), lambda1 = 0.006; end
if nargin < 15 || isempty(lambda2), lambda2 = 0.004; end
if nargin < 16 || isempty(gamma), gamma = 0.012; end
if nargin < 17 || isempty(max_terms_per_equation), max_terms_per_equation = 6; end
if nargin < 18 || isempty(max_quadratic_base_terms), max_quadratic_base_terms = 4; end
if nargin < 19 || isempty(list_anchor_idx)
    list_anchor_idx = ones(size(list_bin_widths));
end
if nargin < 20 || isempty(include_constant_term)
    include_constant_term = false;
end
if nargin < 21 || isempty(include_time_term)
    include_time_term = false;
end
if nargin < 22 || isempty(force_target_peer_terms)
    force_target_peer_terms = false;
end
if nargin < 23 || isempty(max_peer_l3_per_equation)
    max_peer_l3_per_equation = 2;
end

m_interval = plot_end_idx - plot_start_idx + 1;

xobs = [];
mode_labels = {};
mode_levels = [];
Lmax = size(list_modes, 1);

for lev = 1:Lmax
    J = 2^(lev-1);

    for bin = 1:J
        if isempty(list_modes{lev, bin})
            continue;
        end

        eigs_slow = list_w{lev, bin};
        b = list_b{lev, bin};
        bin_width = list_bin_widths(lev, bin);
        t_start = list_t_start(lev, bin);
        t_end = t_start + bin_width - 1;

        if bin_width <= 1 || t_end > m
            continue;
        end

        % Keep only bins overlapping chosen interval
        if t_start > plot_end_idx || t_end < plot_start_idx
            continue;
        end

        % Clip this bin to the selected interval
        local_start = max(t_start, plot_start_idx);
        local_end = min(t_end, plot_end_idx);

        anchor_idx = list_anchor_idx(lev, bin);
        anchor_frame = t_start + anchor_idx - 1;

        local_start = max(local_start, anchor_frame);
        if local_start > local_end
            continue;
        end

        rel_start = local_start - anchor_frame;
        rel_end = local_end - anchor_frame;

        insert_start = local_start - plot_start_idx + 1;
        insert_end = local_end - plot_start_idx + 1;

        for mode_idx = 1:length(eigs_slow)
            a_local = b(mode_idx) * eigs_slow(mode_idx).^(rel_start:rel_end);

            a_full = zeros(m_interval, 1);
            a_full(insert_start:insert_end) = real(a_local(:));

            xobs = [xobs, a_full]; %#ok<AGROW>
            mode_labels{end+1} = sprintf('L%d B%d Mode%d', lev, bin, mode_idx); %#ok<AGROW>
            mode_levels(end+1) = lev; %#ok<AGROW>
        end
    end
end

tobs = (0:m_interval-1)' * dt;

keep_cols = false(1, length(mode_levels));

for lev = input_levels
    lev_cols = find(mode_levels == lev);
    if isempty(lev_cols)
        continue;
    end

    lev_energy = rms(xobs(:, lev_cols), 1);
    [~, ord] = sort(lev_energy, 'descend');
    keep_count = min(top_input_modes_per_level, length(lev_cols));
    keep_cols(lev_cols(ord(1:keep_count))) = true;
end

target_level_cols = find(mode_levels == target_level);
target_energy = rms(xobs(:, target_level_cols), 1);
[~, target_ord] = sort(target_energy, 'descend');
keep_target_count = min(top_target_modes, length(target_level_cols));
keep_cols(target_level_cols(target_ord(1:keep_target_count))) = true;

xobs = xobs(:, keep_cols);
mode_labels = mode_labels(keep_cols);
mode_levels = mode_levels(keep_cols);
target_cols = find(mode_levels == target_level);

weights = [];
polys = 1;
trigs = [];
custom_tags = [];
custom_fcns = {};
phi_class = 1;
max_d = 1;
tau = 10^-16;
tauhat = -1;
K_frac = min(50, floor(m_interval/2));
overlap_frac = 1;
relax_AG = 0;
scale_Theta = 2;
useGLS = 0;

include_quadratic_terms = true;
include_cross_terms = true;
max_active_inputs_per_equation = 4;
relative_term_tol = 0.03;
max_nonlinear_terms_per_equation = 4;

alpha_loss = 0.8;
overlap_frac_ag = 0.8;
mt_ag_fac = [1;0];
pt_ag_fac = [1;0];
run_sindy = 0;
useFD = 2;
smoothing_window = 0;

% After the loop where xobs is constructed:
% Normalize each column of xobs
mu_xobs = mean(xobs, 1);
sigma_xobs = std(xobs, 0, 1);
sigma_xobs(sigma_xobs == 0) = 1;

xobs = (xobs - mu_xobs) ./ sigma_xobs;

[w_linear,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    wsindy_ode_fun_capped(xobs,tobs,weights, ...
    polys,trigs,custom_tags,custom_fcns, ...
    phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG, ...
    scale_Theta,useGLS,lambda1,gamma,alpha_loss, ...
    overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,smoothing_window);

w_second = cell(length(target_cols),1);
labels_second = cell(length(target_cols),1);

n_base = size(xobs,2);

for kk = 1:length(target_cols)

    target_col = target_cols(kk);

    % ---- active input terms from first pass ----
    coef1 = w_linear(:, target_col);
    tol1 = 1e-6 * max(abs(coef1));

    active_linear = find(abs(coef1(1:n_base)) > tol1);

    % Keep a small number of ancestor/input variables as forcing.
    active_inputs = active_linear(ismember(mode_levels(active_linear), input_levels));

    if ~isempty(active_inputs)
        [~, ord] = sort(abs(coef1(active_inputs)), 'descend');
        active_inputs = active_inputs(ord(1:min(max_active_inputs_per_equation, length(active_inputs))));
    end

    active_peers = active_linear(ismember(mode_levels(active_linear), target_level));
    active_peers(active_peers == target_col) = [];

    if force_target_peer_terms
        forced_peers = target_cols(:).';
        forced_peers(forced_peers == target_col) = [];
        active_peers = unique([forced_peers(:); active_peers(:)], 'stable');
    elseif ~isempty(active_peers)
        [~, ord] = sort(abs(coef1(active_peers)), 'descend');
        active_peers = active_peers(ord(1:min(max_peer_l3_per_equation, length(active_peers))));
    end

    if force_target_peer_terms && isfinite(max_peer_l3_per_equation) && length(active_peers) > max_peer_l3_per_equation
        active_peers = active_peers(1:max_peer_l3_per_equation);
    end

    % ---- build second-pass observed state matrix ----
    xobs2 = [];
    labels2 = {};

    if include_constant_term
        xobs2 = [xobs2, ones(size(xobs, 1), 1)]; %#ok<AGROW>
        labels2{end+1} = '1'; %#ok<AGROW>
    end

    if include_time_term
        xobs2 = [xobs2, normalized_time_feature(tobs)]; %#ok<AGROW>
        labels2{end+1} = 't'; %#ok<AGROW>
    end

    % active first-pass linear variables only
    for a = 1:length(active_inputs)
        ii = active_inputs(a);

        xobs2 = [xobs2, xobs(:,ii)]; %#ok<AGROW>
        labels2{end+1} = mode_labels{ii}; %#ok<AGROW>
    end

    % selected peer L3 couplings, kept small for identifiability
    for a = 1:length(active_peers)
        ii = active_peers(a);

        xobs2 = [xobs2, xobs(:,ii)]; %#ok<AGROW>
        labels2{end+1} = mode_labels{ii}; %#ok<AGROW>
    end

    % add target itself
    xobs2 = [xobs2, xobs(:,target_col)];
    labels2{end+1} = mode_labels{target_col};
    target_col2 = length(labels2);

    base_for_quadratics = [active_inputs(:); active_peers(:); target_col];
    base_for_quadratics = unique(base_for_quadratics, 'stable');

    if length(base_for_quadratics) > max_quadratic_base_terms
        base_for_quadratics = base_for_quadratics(1:max_quadratic_base_terms);
    end

    if include_quadratic_terms
        for a = 1:length(base_for_quadratics)
            ii = base_for_quadratics(a);

            xobs2 = [xobs2, xobs(:,ii).^2]; %#ok<AGROW>
            labels2{end+1} = sprintf('(%s)^2', mode_labels{ii}); %#ok<AGROW>
        end
    end

    if include_cross_terms
        for a = 1:length(base_for_quadratics)
            for b = a+1:length(base_for_quadratics)
                ii = base_for_quadratics(a);
                jj = base_for_quadratics(b);

                xobs2 = [xobs2, xobs(:,ii).*xobs(:,jj)]; %#ok<AGROW>
                labels2{end+1} = sprintf('(%s)*(%s)', mode_labels{ii}, mode_labels{jj}); %#ok<AGROW>
            end
        end
    end

    % ---- second pass is still linear in the augmented columns ----
    [w2,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        wsindy_ode_fun_capped(xobs2,tobs,weights, ...
        polys,trigs,custom_tags,custom_fcns, ...
        phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG, ...
        scale_Theta,useGLS,lambda2,gamma,alpha_loss, ...
        overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,smoothing_window);

    coef2 = w2(:, target_col2);

    coef_cap = 25;
    if any(abs(coef2) > coef_cap)
        fprintf('Clipping %s coefficients above %.1f in reduced model.\n', ...
            mode_labels{target_col}, coef_cap);
        coef2 = max(min(coef2, coef_cap), -coef_cap);
    end

    self_idx = find(strcmp(labels2, mode_labels{target_col}), 1);
    %{
    if ~isempty(self_idx) && coef2(self_idx) > 0
        fprintf('Removing positive self-feedback for %s: %.4e -> 0\n', ...
            mode_labels{target_col}, coef2(self_idx));
        coef2(self_idx) = 0;
    end
    %}
    self_growth_cap = 0.02;

    if ~isempty(self_idx) && coef2(self_idx) > self_growth_cap
        fprintf('Capping positive self-feedback for %s: %.4e -> %.4e\n', ...
            mode_labels{target_col}, coef2(self_idx), self_growth_cap);
        coef2(self_idx) = self_growth_cap;
    end

    coef2 = prune_equation_terms(coef2, labels2, mode_labels{target_col}, ...
        relative_term_tol, max_terms_per_equation, max_nonlinear_terms_per_equation);

    w_second{kk} = coef2;
    labels_second{kk} = labels2;
end

end

function t_feature = normalized_time_feature(tobs)
%NORMALIZED_TIME_FEATURE Scale local time to [-1, 1].

t0 = min(tobs);
t1 = max(tobs);
if t1 <= t0
    t_feature = zeros(size(tobs));
else
    t_feature = 2 * (tobs - t0) / (t1 - t0) - 1;
end

end





