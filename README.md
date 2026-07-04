# mrDMD + WSINDy CFD Experiments

This folder organizes the original exploratory MATLAB scripts into a small project.

## Main entry points

- `main/run_single_experiment.m` runs one configured mrDMD + WSINDy experiment.
- `main/run_parameter_sweep.m` runs the WSINDy tuning sweep.
- `scripts/data_prep/create_cfd_mat.py` creates the cropped/downsampled `.mat`
  file from extracted CFD `.txt` snapshots.

## Folder layout

- `config/` contains baseline settings and sweep grids.
- `config/datasets/` contains one runnable config per dataset.
- `src/data/` contains CFD dataset loading and validation.
- `src/mrdmd/` contains DMD, mrDMD, reconstruction, and modal utility code.
- `src/wsindy/` contains the mrDMD-to-WSINDy workflow and evaluation code.
- `src/visualization/` contains plotting helpers.
- `results/` is for generated sweep tables, figures, and fitted models.
- `data/raw/` is for raw or copied CFD text snapshots if you want them in-project.
- `data/processed/` is for generated `.mat` files.
- `tests/` is reserved for small validation scripts.

## Common plot commands

## Choosing A Dataset

The pipeline expects every dataset to provide `X` as `states x time`.

For the current CFD velocity dataset:

```matlab
clear cfg
startup
cfg = cfd_config();
cfg.data.plot_field = 'speed';  % use 'u', 'v', or 'speed'
run('main/run_single_experiment.m')
```

For the ducks scalar dataset:

```matlab
clear cfg
startup
cfg = ducks_config();
run('main/run_single_experiment.m')
```

For the earthquake scalar dataset:

```matlab
clear cfg
startup
cfg = earthquake_config();
run('main/run_single_experiment.m')
```

Scalar `.mat` files should contain:

```matlab
X
```

and preferably either:

```matlab
x, y
```

or:

```matlab
image_shape   % [rows cols]
```

If no coordinates are supplied, the loader falls back to index-based plotting.

The dataset-specific config files live here:

```text
config/datasets/cfd_config.m
config/datasets/ducks_config.m
config/datasets/earthquake_config.m
```

Update the `cfg.frames.*` values inside the ducks and earthquake config files
once those datasets are copied into `data/processed/` and you know the snapshot
counts.

## Result Saving

Result behavior is controlled in `config/default_experiment_config.m`:

```matlab
cfg.results.save_single_model = true;
cfg.results.save_sweep_results = true;

cfg.results.clear_models_on_single_run = false;
cfg.results.clear_sweeps_on_sweep_run = false;
```

Set the `clear_*` switches to `true` only when you intentionally want to delete old generated result files in that result folder before a run.

After `main/run_single_experiment.m`:

```matlab
plot_l3_error(details, cfg)
plot_full_error(details, cfg)
plot_l3_frame(details, data, 1)
plot_full_frame(details, data, 1)
plot_experiment_summary(details, data, cfg, 1)
```

Animation commands:

```matlab
animate_l3_wsindy_comparison(details, data, cfg)
animate_full_wsindy_comparison(details, data, cfg)
animate_mrdmd_mode_groups(mrdmd, data)
animate_experiment_summary(details, data, cfg, mrdmd)
```

After `main/run_parameter_sweep.m`:

```matlab
plot_sweep_results(result_table)
```

## Notes

The code still expects the WSINDy dependency functions used by your original
`wsindy_ode_fun` implementation to be on the MATLAB path, including functions
such as `get_tags`, `build_theta`, `sparsifyDynamics`, `phi_int_weights`,
`findcorners`, and related WSINDy utilities.

By default, the MATLAB scripts load:

`data/processed/cylinder_2D_uv_ONE_MIDDLE_Z_cropped_downsampled.mat`

To generate it from your extracted text snapshots, run:

```powershell
python scripts/data_prep/create_cfd_mat.py --data-dir "C:\Users\edric\Desktop\2d_cfd_data"
```

Or update `cfg.data.mat_file` in `config/default_experiment_config.m`.
