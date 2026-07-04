function data = validate_spatiotemporal_dataset(raw, cfg)
%VALIDATE_SPATIOTEMPORAL_DATASET Validate datasets with spatial states over time.
%
% Supported layouts:
%   uv     : X = [u; v], states x time, size(X,1) = 2*npoints.
%   scalar : X = scalar spatial state, states x time, size(X,1) = npoints.
%   auto   : infer uv when possible, otherwise scalar.

if ~isfield(raw, 'X')
    error('The dataset should include X as states x time.');
end

X = raw.X;
layout = lower(cfg.data.state_layout);

if strcmp(layout, 'auto')
    if isfield(raw, 'npoints') && size(X,1) == 2*double(raw.npoints)
        layout = 'uv';
    else
        layout = 'scalar';
    end
end

switch layout
    case 'uv'
        [npoints, x, y] = validate_uv_layout(raw, X);
        state_components = 2;
        component_names = {'u', 'v'};

    case 'scalar'
        [npoints, x, y] = validate_scalar_layout(raw, X);
        state_components = 1;
        component_names = {'scalar'};

    otherwise
        error('Unknown cfg.data.state_layout "%s". Use ''uv'', ''scalar'', or ''auto''.', cfg.data.state_layout);
end

if isfield(raw, 'dt')
    dt = raw.dt;
else
    dt = 1; % Use snapshot index as time if the MAT file does not include physical dt.
end

data.X = X;
data.x = x;
data.y = y;
data.dt = dt;
data.npoints = npoints;
data.state_layout = layout;
data.plot_field = lower(cfg.data.plot_field);
data.state_components = state_components;
data.component_names = component_names;

if isfield(raw, 'z')
    data.z = raw.z(:);
end
if isfield(raw, 'tvals')
    data.tvals = raw.tvals(:);
end
if isfield(raw, 'image_shape')
    data.image_shape = raw.image_shape;
end

end

function [npoints, x, y] = validate_uv_layout(raw, X)
if isfield(raw, 'npoints')
    npoints = double(raw.npoints);
else
    if mod(size(X,1), 2) ~= 0
        error('Expected X to be stacked as [u; v], but size(X,1) is odd.');
    end
    npoints = size(X,1) / 2;
end

if size(X,1) ~= 2*npoints
    error('Expected size(X,1) to equal 2*npoints for stacked [u; v] snapshots.');
end

if ~isfield(raw, 'x') || ~isfield(raw, 'y')
    error('The [u; v] dataset should include point coordinates x and y.');
end

x = raw.x(:);
y = raw.y(:);
if numel(x) ~= npoints || numel(y) ~= npoints
    error('Expected x and y to each contain npoints entries.');
end
end

function [npoints, x, y] = validate_scalar_layout(raw, X)
if isfield(raw, 'npoints')
    npoints = double(raw.npoints);
else
    npoints = size(X,1);
end

if size(X,1) ~= npoints
    error('Expected size(X,1) to equal npoints for scalar spatial snapshots.');
end

if isfield(raw, 'x') && isfield(raw, 'y')
    x = raw.x(:);
    y = raw.y(:);
elseif isfield(raw, 'image_shape')
    image_shape = double(raw.image_shape(:));
    if numel(image_shape) ~= 2 || prod(image_shape) ~= npoints
        error('image_shape must be [rows cols] and prod(image_shape) must equal npoints.');
    end
    [xx, yy] = meshgrid(1:image_shape(2), 1:image_shape(1));
    x = xx(:);
    y = yy(:);
else
    % Fallback for scalar datasets that have spatial states but no coordinates.
    % This preserves the analysis path; plots become index-based strips.
    x = (1:npoints).';
    y = zeros(npoints, 1);
end

if numel(x) ~= npoints || numel(y) ~= npoints
    error('Expected x and y to each contain npoints entries.');
end
end

