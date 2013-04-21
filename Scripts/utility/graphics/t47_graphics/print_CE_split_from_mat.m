function print_CE_split_from_mat(mat_file, source, space)

%PRINT_CE_SPLIT(mat_file, source, [h w sp ln tl])
%
%by T47, Feb 2013
%
%Prints the mutate-and-map electrophoregram in h*w splitted figures, with X-axis
%   labelled with mutants' name, Y-axis labelled with each nucleotides' position.
%
%Arguments
%=========
%   mat_file        Required    Provides the .mat file location.
%   source          Required    Provides the source of certain arguments (seqpos,
%                                   sequence, offset).
%                               0 equals from .mat file directly, 1 equals from
%                                   d_rdat inside .mat file.
%   [h w sp ln tl]  Optional    Provides the layout format that will split into h*w 
%                                   pages with vertical blank edge sp. 
%                               h denotes the number of pages vertically;
%                               w denotes the number of pages horizontally;
%                               ln denotes whether to make lines every 12 lanes or not;
%                               tl denotes whether title is added to figure.
%                                   0 equals false; 1 equals true.
%                               Default [2 2 50 1 1]: print in 2*2, spacer 50, with
%                                   lines, with title.
%
% e.g. PRINT_CE_SPLIT(mat_file, 0, []);
%                                   
%
%Notes
%=====
%Print all the figures on US Letter, with "landscape fit":
%   [left -0.25; top 0.25; width 11.50; height 8.50]
%   For bottom figures, use top 0 to print X-axis fully; for top figures, use top 0.50
%       to print title fully.
%Spaces on each border are included for easy splicing.
%All current opened figures will be lost. Save before run.
%

if source == 0
    load(mat_file, 'd_align', 'd_rdat', 'seqpos', 'xsel', 'sequence', 'offset');
else
    load(mat_file, 'd_align', 'd_rdat', 'xsel');
    seqpos = d_rdat.seqpos;
    sequence = d_rdat.sequence;
    offset = d_rdat.offset;
end;

print_CE_split(d_align, d_rdat, seqpos, xsel, sequence, offset, space);

end
