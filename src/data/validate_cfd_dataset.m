function data = validate_cfd_dataset(raw)
%VALIDATE_CFD_DATASET Validate expected CFD MAT-file fields.

cfg.data.state_layout = 'uv';
cfg.data.plot_field = 'speed';
data = validate_spatiotemporal_dataset(raw, cfg);

end
