function [ data_all, filenames_all, data_init, data_length ] = ...
    read_abi_dirs_GUI( filepath, dirnames );
% Read data
% From original code, we add two lines in order to count number of capillaries

data = {};
filenames_all = {};
count = 0;
for j = 1:length( dirnames )
  % This script calls read_abi.m which has the actual file format.
  fprintf( 1, 'Reading in...%s\n',dirnames{j} );
      [data,filenames] = plot_ABI_runs_GUI( [filepath, dirnames{j}], 1 ); 
  
  % Count number of capillaries 
  data_init(j) = count + 1;           % acl
  data_length(j) = length(filenames); % acl
  
  for k = 1:length( filenames )
    count = count + 1;
    data_all{ count }  = data{ k };
      filenames_all{ count }  = filenames{ k };
  end
end


clear count;
clear j; clear k;
