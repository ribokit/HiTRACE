function []=varna_fig(filename,sequence,structure,DATA,colorscheme,offset,special_base_pairs,special_colors, bpp_values, bpp_anchor_bases)
% VARNA_FIG: Create html file with secondary structure -- double click to get VARNA visualizer in a web browser 
%
%  varna_fig(filename,sequence,structure,DATA,colorscheme,offset,special_base_pairs,special_colors, bpp_values, bpp_anchor_bases)
%
% filename  = output filename [e.g., 'my_rna_in_varna.html']
% sequence  = RNA sequence
% structure = structure in dot/bracket notation. length should match sequence.
% 
% optional input parameters:
% DATA               = data for coloring residues. will only use values between 0 and 1, 
%                        so normalize ahead of time. Give [] if no data.
% colorscheme        = blue/white/red (0), another blue/white/red (1), white/orange/red (2)
% special_base_pairs = extra pairs of residues that should have connecting lines.
% special_colors     = colors of connecting lines for special_base_pairs.
% bpp_values         = numbers that will show up on helices as percentages
% bpp_anchor_bases   = where those bpp_values will show up.
%
% (C) R. Das 2011
% (C) C.C. VanLang, P. Cordero 2010

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'colorscheme','var' ); colorscheme = 1; end;

if ~isempty( DATA )
  %% DATA Prep
  %%%Z=hSHAPE(DATA);
  Z=DATA;
  
  Z = max( Z, 0 );
  Z = min( Z, 2.0 );
  
  %Z=100*Z;

  graypoints = find( DATA == -999 | isnan( DATA ) );
  Z( graypoints ) = -0.01;

end



%%
%filename='VARNA/test.html';
fid=fopen(filename,'w');

vers=1;


fprintf(fid,'%s\n','<HTML><HEAD><META http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></HEAD><BODY>');
fprintf(fid,'%s\n','<APPLET code="VARNA.class" codebase="http://rmdb.stanford.edu/site_media/bin" archive="VARNA.jar" width="1200" height="1200">');
fprintf(fid,'%s%s%s\n','<PARAM name="sequenceDBN"  value="',sequence,'"/>');
fprintf(fid,'%s%s%s\n','<PARAM name="structureDBN" value="',structure,'"/>');
fprintf(fid,'%s\n','<PARAM name="algorithm" value="radiate" />');

if ~isempty( DATA )
  fprintf(fid,'%s','<param name="colorMap" value="');
  for i=1:length(Z)
    fprintf(fid,' %6.3f%s',Z(i),',');
  end
  fprintf(fid,'%s\n','"/>');
end

if exist( 'Z','var' )
  switch colorscheme
   case 0 % previous default
    fprintf(fid,'%s\n','<param name="colorMapStyle" value="0:#0000FF;10:#0000FF;40:#FFFFFF;60:#FFFFFF;90:#FF0000;100:#FF0000" />'); 
   case 1 % new default
    fprintf(fid,'%s\n','<param name="colorMapStyle" value="0:#0000FF;0.5:#FFFFFF;1:#FF0000" />'); 
   case 2 % white orange to red
	  % slight pain because of VARNA rescaling:
	  if ( sum( Z < 0 ) > 0 )
	    fprintf(fid,'%s\n','<param name="colorMapStyle" value="-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000" />'); 
	  else
	    fprintf(fid,'%s\n','<param name="colorMapStyle" value="0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000" />'); 
	  end
   case 3 % blue to white to yellow
	  % slight pain because of VARNA rescaling:
	  if ( sum( Z < 0 ) > 0 )
	    fprintf(fid,'%s\n','<param name="colorMapStyle" value="-0.001:#B0B0B0,0:#0000FF;0.5:#FFFFFF;1:#FFFF00" />'); 
	  else
	    fprintf(fid,'%s\n','<param name="colorMapStyle" value="0:#FFFFFF;,0.5:#2222FF;1:#0000FF" />'); 
	  end
  end
end

fprintf( fid, '<param name="bpStyle" value="simple" />\n' );
fprintf( fid, '<param name="baseInner" value="#FFFFFF" />\n' );
fprintf( fid, '<param name="baseOutline" value="#FFFFFF" />\n' );
fprintf( fid, '<param name="bp" value="#000000" />\n' );

if exist( 'special_base_pairs','var' )
  if length( special_base_pairs ) ~= length( special_colors );  fprintf( 'Must specify a special_color for each special_base_pair set\n'); end;

  fprintf(fid,'<param name="auxBPs" value="' ); 
  
  for q = 1:length( special_base_pairs )
    special_base_pair_set = special_base_pairs{q};
    hex_color =  convert_rgb_to_hexadecimal( special_colors{q} );
    for k = 1:size( special_base_pair_set, 1);
      fprintf( fid, '(%d,%d):thickness=3,color=#%6s;', special_base_pair_set(k,1), special_base_pair_set(k,2), hex_color );
    end
  end
  fprintf( fid, '" />\n' );
    
end

if ( exist( 'offset','var' ) || exist( 'bpp_values','var' ) )
  
  if exist( 'offset','var' )
    fprintf( fid, '<param name="baseNum" value="#FFFFFF" />\n' );
    fprintf( fid, '<param name="periodNum" value="1000" />\n' );
  end
  
  fprintf( fid, '<param name="annotations" value="' );
  if exist( 'offset','var' )
    PERIOD = 10;
    for i = 1:length( sequence )
      if ( mod( i+offset, PERIOD) == 0 )
	fprintf( fid, '%d:type=B,anchor=%d,color=#000000,size=8;', i+offset, i );
      end
    end
  end

  if exist( 'bpp_values','var' ) && ~isempty( bpp_values )
    if length( bpp_values ) ~= length( bpp_anchor_bases );  fprintf( 'Must specify a bpp_anchor_base for each bpp_value \n'); end;
    for i = 1:length( bpp_values )
      %fprintf( fid, '%3.0f%%:type=L,anchor=%d,color=#303030,size=9;', 100*bpp_values(i), bpp_anchor_bases(i) );
      fprintf( fid, '%3.0f%%:type=L,anchor=%d,color=#FF3030,size=9;', 100*bpp_values(i), bpp_anchor_bases(i) );
    end
  end
   

  fprintf( fid, '">\n' );

end




fprintf(fid,'%s\n','<PARAM name="flat" value="true" />');
fprintf(fid,'%s\n','</applet></BODY></HTML>');

fclose( fid );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hex_string = convert_rgb_to_hexadecimal( rgb )

hex_string = '';
for i = 1:3
  hex_string = [ hex_string, pad_with_zero( dec2hex( floor(rgb(i)*255) ) ) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = pad_with_zero( s )
if length( s) == 1; 
  s = ['0',s];
end