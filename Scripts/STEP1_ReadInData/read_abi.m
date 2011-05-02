function data = read_abi( filename );
% READ_ABI
%
% data = read_abi( filename )
%    Returns four data traces from ABI 3100 (or generally
%     .ab1 format) files. 
%    NB: There's a lot more information in each file, e.g.,
%     voltages, wattage, number of dyes, etc. More could be
%     extracted.
%
%    The file format is base on http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
%
%  Rhiju Das, Nov. 2008.
%

%filename = ['Das1_Run_3100_2008-11-03_455/D03_X20-Cy3-ddTTP_07.ab1'];
%filename = ['Das1_Run_3100_2008-11-03_455/D04_X20-Cy3-ddTTP-dil_08.ab1'];
%filename = ['Das1_Run_3100_2008-11-03_455/H03_X20-Mix-TAMRA-ddGTP-Cy3-ddTTP_15.ab1'];

[fid,message] = fopen( filename, 'r', 'b'); 
dir.abi = fread( fid, 4, '*char' );
dir.version = fread( fid, 1, 'uint16');
dir.name = fread( fid, 4, '*char' );
dir.number = fread( fid, 1, 'int32' );
dir.elementtype = fread( fid, 1, 'int16' );
dir.elementsize = fread( fid, 1, 'int16' );
dir.numelements = fread( fid, 1, 'int32' );
dir.datasize = fread( fid, 1, 'int32' );
dir.dataoffset = fread( fid, 1, 'int32' );
dir.datahandle = fread( fid, 1, 'int32' );
dir.dummy =  fread( fid, 47, 'int16' ) ;

status = fseek( fid, dir.dataoffset, -1 );

name = [];
for j = 1: dir.numelements 
  name(:,j) = fread( fid, 4, '*char')' ;
  number(j) = fread( fid, 1, 'int32');
  elementtype(j) = fread( fid, 1, 'int16');
  elementsize(j) = fread( fid, 1, 'int16');
  numelements(j) = fread( fid, 1, 'int32');
  datasize(j) = fread( fid, 1, 'int32');
  dataoffset(j) = fread( fid, 1, 'int32');
  datahandle(j) = fread( fid, 1, 'int32');
end
name = char( name )';

count = 0;
data_terms = {};
for j = 1:dir.numelements
  if ( strfind( 'DATA', name(j,:)) )
    count = count + 1;
    status = fseek( fid, dataoffset(j), -1 );
    data_terms{ count } = fread( fid, numelements(j), 'int16' )';
  end
end

data = [];
for j = 1:4
  %data(:,j) = smooth(data_terms{j});
  data(:,j) = data_terms{j};
  data(:,j) = baseline_subtract( data(:,j) );
end


fclose( fid );
