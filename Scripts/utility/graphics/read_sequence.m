function sequence = read_sequence( filename );
%  sequence = read_sequence( filename );

if nargin == 0;  help( mfilename ); return; end;

fid = fopen( filename );
sequence = fgetl(fid);
fprintf(1,'%s\n',sequence );
fclose(fid);
