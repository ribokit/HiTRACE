function sequence = read_sequence( filename );
% READ_SEQUENCE
%
%  sequence = read_sequence( filename );
%
% simple read in of one line from a text file.
%
fid = fopen( filename );
sequence = fgetl(fid);
fprintf(1,'%s\n',sequence );
fclose(fid);
