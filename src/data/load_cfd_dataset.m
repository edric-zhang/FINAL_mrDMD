function data = load_cfd_dataset(cfg)
%LOAD_CFD_DATASET Compatibility wrapper for older scripts.

cfg.data.state_layout = 'uv';
if ~isfield(cfg.data, 'plot_field')
    cfg.data.plot_field = 'speed';
end
data = load_spatiotemporal_dataset(cfg);

end
