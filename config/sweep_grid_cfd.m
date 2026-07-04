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
        % lambda1  lambda2  gamma   top_inputs  top_targets  max_terms  quad_base

        % Top targets = 3, conservative
        0.30      0.06     0.001     3      3      3      2
        0.40      0.08     0.001     3      3      4      3
        0.50      0.10     0.001     3      3      5      4

        % Top targets = 3, current neighborhood
        0.30      0.06     0.001     5      3      5      4
        0.40      0.08     0.001     5      3      5      4
        0.50      0.10     0.001     5      3      6      5

        % Top targets = 3, more inputs / richer
        0.20      0.04     0.001     8      3      6      5
        0.40      0.08     0.001     8      3      8      6
        0.80      0.16     0.001     8      3      8      6

        % Top targets = 3, gamma sensitivity
        0.40      0.08     0.005     5      3      5      4
        0.40      0.08     0.010     5      3      5      4
        0.40      0.08     0.050     5      3      5      4

        % Top targets = 4, conservative
        0.30      0.06     0.001     3      4      3      2
        0.40      0.08     0.001     3      4      4      3
        0.50      0.10     0.001     3      4      5      4

        % Top targets = 4, current neighborhood
        0.30      0.06     0.001     5      4      5      4
        0.40      0.08     0.001     5      4      5      4
        0.50      0.10     0.001     5      4      6      5

        % Top targets = 4, more inputs / richer
        0.20      0.04     0.001     8      4      6      5
        0.40      0.08     0.001     8      4      8      6
        0.80      0.16     0.001     8      4      8      6

        % Top targets = 4, gamma sensitivity
        0.40      0.08     0.005     5      4      5      4
        0.40      0.08     0.010     5      4      5      4
        0.40      0.08     0.050     5      4      5      4
        ]; 
    
    %{
    param_grid = [
    % lambda1 lambda2 gamma top_inputs top_targets max_terms quad_base
    0.4 0.08 0.001  5  1  5  4
    0.4 0.08 0.001  5  2  5  4
    0.4 0.08 0.001  5  3  5  4
    0.4 0.08 0.001  5  4  5  4
    0.4 0.08 0.001  5  6  5  4
    0.4 0.08 0.001  5  8  5  4
    ];
    %}
    
    %{



    clear; clc; close all;
    startup;
    cfg = earthquake_config();
    cfg.results.save_sweep_results = false;
    cfg.results.clear_sweeps_on_sweep_run = false;
    run('main/run_parameter_sweep.m');

    clear; clc; close all;
    startup;
    cfg = earthquake_config();
    cfg.results.save_sweep_results = false;
    cfg.results.clear_sweeps_on_sweep_run = false;
    run('main/run_single_experiment.m');

    %}

end

