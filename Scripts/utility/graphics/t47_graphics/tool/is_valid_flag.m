function flag_out = is_valid_flag(flag_in, flag_default)

%
% flag_out = IS_VALID_FLAG(flag_in, flag_default);
%
% Checks input flag array and supplement the missing ones from
%  flag_default. If input has more elements than default, no changes will
%  be made.
%
% by T47, May 2013.
%

flag_out = flag_in;
if length(flag_in) < length(flag_default);
    flag_default(1:length(flag_in)) = flag_in;
    flag_out = flag_default;
end;

if length(flag_in) > length(flag_default);
    fprintf('WARNING: input has more elements than default.\n');
end;

