function mode_groups = default_mode_groups(data)
%DEFAULT_MODE_GROUPS Default grouped mrDMD modes from the original script.

if nargin >= 1 && isfield(data, 'state_layout') && strcmpi(data.state_layout, 'scalar')
    mode_groups = {
        [3 4 2],                  'scalar', '3-4-2: Scalar contribution'
        [2 2 11],                 'scalar', '2-2-11: Scalar contribution'
        [3 4 2; 2 2 11],          'scalar', '3-4-2 & 2-2-11: Combined scalar contribution'
    };
else
    mode_groups = {
        %[1 1 1; 1 1 2],          'u', 'L1 slow/global envelope';
        %[1 1 7],                 'u', '1-1-7: Primary Von Karman Vortex Street: SHREDDING';
        [3 4 2],                  'u', '3-4-2: Medium-size, wavelike ripplings (u)'
        [2 2 11],                 'u', '2-2-11: Broad, elongated shape (u)'
        [2 2 11],                 'v', '2-2-11: Broad, elongated shape (v)'
        [3 4 2; 2 2 11],          'u', '3-4-2 & 2-2-11: Both: Shred --> Shear-Layer Transition (u)'
        %[2 2 12; 3 4 4],         'u', '2-2-12 & 3-4-4: KH Shear Layer Instabilities';
    };
end

end
