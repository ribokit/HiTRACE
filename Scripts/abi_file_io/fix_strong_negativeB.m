function d = fix_strong_negativeB(data,range)

% The 350-ROX ladder from GeneScan contains a strong band that migrates at 
% a smaller molecular weight than the 35-nucleotide marker. This band is
% almost always saturating and cannot be corrected by the leakage
% correction in quick_look, leading to problems with the following smooth
% baseline subtraction. Thus, use interpolation to flatten the strong
% negative peak before the smooth baseline subtraction step in quick_look.
%
% Clarence Cheng, 2013

xi = range;
if range(1) == 0 | isempty( range ); d = data; return; end;
  
frame_pos = [xi(1); xi(length(xi))];

% flatten
frame_val = [data(frame_pos(1),2); data(frame_pos(2),2)];     %define intensity values at timepoints framing the to-be-flattened region
yi = interp1(frame_pos,frame_val,xi);      %return interpolated values between framing timepoints of the negative peak into yi
data2 = data;
for i = 1:length(xi)
  data2(xi(i),2) = yi(i);     %replace values in data with interpolated values
end

d = data2;

%figure(6); plot(data(:,2),'b'); hold on; plot(data2(:,2),'Color',[0 0.5 0]);     %plot uninterpolated and interpolated data