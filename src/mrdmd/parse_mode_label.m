function [lev, bin, mode_idx] = parse_mode_label(label)
%PARSE_MODE_LABEL Turns something like 'L3 B4 Mode14' into numeric indices.

nums = sscanf(label, 'L%d B%d Mode%d');

lev = nums(1);
bin = nums(2);
mode_idx = nums(3);

end
