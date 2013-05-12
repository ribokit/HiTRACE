function works = parallelization_exists();
%  PARALLELIZATION_EXISTS:  helper function to handle matlabpool parallelization.
%
% works = parallelization_exists();
%
%

works = 0;

if exist( 'matlabpool' )  
  works = 1;
  try
    if matlabpool( 'size' ) == 0 ;   
      fprintf( 'Starting up MATLAB parallelization toolbox...\n' );
      res = findResource; 
      num_daemons = round(res.ClusterSize);
      fprintf( 'Going to try %d daemons... if you want more use Parallel/Manage Configurations... to create a new local configuration\n', num_daemons );      
      matlabpool( num_daemons ); 
    end
  catch me
    fprintf( 'WARNING! NOT RUNNING PARALLELIZATION TOOLBOX!\n');
    fprintf( 'Check out: http://www.mathworks.com/support/bugreports/919688\n');
    works = 0;
  end
end



