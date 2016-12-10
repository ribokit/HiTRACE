function flag = check_graphic_interface()
% return 1 if 'old'(<= R2014a) graphic interface is available
% return 0 if not (>= R2014b)

[year, sub_release] = get_matlab_version();
if year > 2014 || (year == 2014 && strcmp(sub_release, 'b'));
    flag = 0;
else
    flag = 1;
end;
