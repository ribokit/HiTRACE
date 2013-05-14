function seqpos = get_seqpos( sequence, offset, xsel )

%seqpos = length( sequence ) + offset + 1 - [1:length(xsel)];
seqpos = length( sequence ) + offset - length(xsel) + [1:length(xsel)];