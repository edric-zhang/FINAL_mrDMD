function param_grid = sweep_grid_cfd()
%SWEEP_GRID_CFD Candidate WSINDy structures/parameters for CFD tuning.
%
% Columns:
% lambda1, lambda2, gamma,
% top_input_modes_per_level, top_target_modes,
% max_terms_per_equation, max_quadratic_base_terms
%
% CFD note:
% Your previous full-field error was already low, but L3 percent error was huge.
% So this sweep varies structure more aggressively than lambda alone.

param_grid = [
    0.004 0.002 0.010   5 10   5 4
    0.006 0.004 0.010   5 10   5 4
    0.008 0.006 0.010   5 10   5 4

    0.004 0.002 0.011   5 10   5 4
    0.006 0.004 0.011   5 10   5 4
    0.008 0.006 0.011   5 10   5 4

    0.004 0.002 0.012   5 10   5 4
    0.006 0.004 0.012   5 10   5 4
    0.008 0.006 0.012   5 10   5 4

    0.004 0.002 0.013   5 10   5 4
    0.006 0.004 0.013   5 10   5 4
    0.008 0.006 0.013   5 10   5 4

    0.004 0.002 0.014   5 10   5 4
    0.006 0.004 0.014   5 10   5 4
    0.008 0.006 0.014   5 10   5 4

    % Check 9 and 11 targets around same gamma
    0.006 0.004 0.011   5 9    5 4
    0.006 0.004 0.012   5 9    5 4
    0.006 0.004 0.013   5 9    5 4

    0.006 0.004 0.011   5 11   5 4
    0.006 0.004 0.012   5 11   5 4
    0.006 0.004 0.013   5 11   5 4

    % Check whether 6 terms helps target 10
    0.006 0.004 0.011   5 10   6 4
    0.006 0.004 0.012   5 10   6 4
    0.006 0.004 0.013   5 10   6 4
];

end

