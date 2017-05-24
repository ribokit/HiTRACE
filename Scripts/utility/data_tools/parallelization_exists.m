function works = parallelization_exists();
%  PARALLELIZATION_EXISTS:  helper function to handle matlabpool parallelization.
%
% works = parallelization_exists();
% 
% makes use of parpool (or on early versions of MATLAB, matlabpool)
% 
% (C) Das lab 2012-2017

works = 0;
if exist( 'parpool' )
    works = 1;
    try
        if isempty(gcp('nocreate'))
            fprintf( 'Starting up MATLAB parallelization toolbox...\n' );
            mycluster = parcluster();
            num_daemons = mycluster.NumWorkers;
            fprintf( 'Going to try %d daemons... if you want more use Parallel/Manage Configurations... to create a new local configuration\n', num_daemons );
            parpool( num_daemons );
        end
    catch me
        fprintf( 'WARNING! NOT RUNNING PARALLELIZATION TOOLBOX!\n');
        fprintf( 'Check out: http://www.mathworks.com/support/bugreports/919688\n');
        works = 0;
    end
    currentpool = gcp;
    fprintf( 'Current worker pool has %d daemons\n', currentpool.NumWorkers );
elseif exist( 'matlabpool' )
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



