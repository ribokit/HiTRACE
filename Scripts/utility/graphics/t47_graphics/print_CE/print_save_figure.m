function print_save_figure(fig, file_name, dir_name, flg)

%
% PRINT_SAVE_FIGURE(fig, file_name, dir_name, flg);
%
% Output figure to .eps prints and MATLAB figure files.
%
% Input
% =====
%   fig         Required        Provides the figure handle, e.g. figure(2)
%                                   or gcf.
%   file_name   Required        Provides the output file name. Format in
%                                   string. Extension will automatically
%                                   append.
%   dir_name    Optional        Provides the output folder name. Format in
%                                   string. It will the create folder if 
%                                   does not exist. Default is 'Figures'.
%   flag        Optional        Provides whether to output to .fig files.
%                                   0 for NO, 1 for YES. Default is 0.
%
%
% by CYC, T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

file_name = strrep(strrep(file_name, '/', ''), ' ', '');

current_dir = pwd;

if ~exist('dir_name', 'var') || isempty(dir_name); dir_name = 'Figures'; end;
dir_name = strcat(current_dir, '/', dir_name);
if ~exist(dir_name, 'dir'); mkdir(dir_name); end;

full_path = [dir_name, '/', file_name,'.eps'];
print(fig, '-depsc2', '-loose', '-r300', full_path);
fprintf( ['Created: ', full_path, '\n'] );

if ~exist('flg','var') || isempty(flg); flg = 0 ; end;
if flg == 1;
    tag = [dir_name, '/', file_name,'.fig'];
    hgsave(fig, tag);
end;

clear dir_name;
