function [ sequence_reversed ] = sequence_reversed( sequence )
%This function reverses an input array

sequence_reversed = [];
for i = 1:length(sequence)
    sequence_reversed(i) = sequence(length(sequence) + 1 - i);
end;

return;

