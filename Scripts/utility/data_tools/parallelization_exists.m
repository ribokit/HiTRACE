function works = parallelization_exists();
% 
% works = parallelization_exists();
%
works = 0;

if exist( 'matlabpool' )  
  works = 1;
  try
    if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end
  catch me
    fprintf( 'WARNING! NOT RUNNING PARALLELIZATION TOOLBOX!\n');
    fprintf( 'Check out: http://www.mathworks.com/support/bugreports/919688\n');
    works = 0;
  end
end



