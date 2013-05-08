function out = is_valid_boolean (in)

out = in;
if (in ~= 0) && (in ~= 1);
    out = 0;
    fprintf('WARNING: Invalid boolean flag input.\n');
end;