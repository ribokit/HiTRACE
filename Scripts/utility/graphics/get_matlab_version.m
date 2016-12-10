function [year, sub_release] = get_matlab_version()

ver_str = version('-release');
year = str2num(ver_str(1:end-1));
sub_release = ver_str(end);
