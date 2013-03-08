function [ sequence_reversed ] = untitled2( sequence )
%This function reverses an input array


sequence_reversed = [];
for i = 1:length(sequence)
    sequence_reversed(i) = sequence(length(sequence) + 1 - i);
end;
end

