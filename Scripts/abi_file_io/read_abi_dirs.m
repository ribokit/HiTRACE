function [ data_all, filenames_all, data_set_starts, data_length ] = ...
    read_abi_dirs( filepath, dirnames, dye_names_full, PLOT_STUFF );
% READ_ABI_DIRS:  wrapper around plot_abi_runs -- read in .abi data from multiple directories.
%
% INPUTS:
%  filepath   = path to directories with ABI files.
%  dirnames   = directory with ABI files [.ab1 or .fsa format]
%  dye_names_full = [optional] names of dyes in each color channel. default = {}, which means no leakage correction. 
%                   Can also specify filename with leakage matrix.
%  PLOT_STUFF = [optional, ignore for now] default: 1.
%
% OUTPUTS:
%  data_all        = cell containing data matrices for all .ab1/.fsa files. One matrix for each color channel.
%  filenames_all   = which .ab1/.fsa files were read in.
%  data_set_starts = index of the start of each data subset [corresponding to a different directory]
%  data_length     = number of files read in from each data subset
%
% (C) R. Das 2008-2011, 2013
%

if nargin == 0;  help( mfilename ); return; end;

data = {};
data_all = {};
data_set_starts = [];
data_length = [];
filenames_all = {};
count = 0;

if ~exist( 'dye_names_full' ) dye_names_full = {}; end;

for j = 1:length( dirnames )
  % This script calls read_abi.m which has the actual file format.
  fprintf( 1, 'Reading in:  %s\n',dirnames{j} ); 
  
  [data,filenames] = plot_ABI_runs( [filepath, dirnames{j}], dye_names_full, PLOT_STUFF ); 

  % Count number of capillaries 
  data_set_starts(j) = count + 1;           % acl
  data_length(j) = length(filenames); % acl
  
  for k = 1:length( filenames )
    count = count + 1;
    data_all{ count }  = data{ k };
    filenames_all{ count }  = filenames{ k };
  end
end


clear count;
clear j; clear k;
