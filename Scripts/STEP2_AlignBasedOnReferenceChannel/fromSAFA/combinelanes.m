function profiles_combine=combinelanes(profiles_align,numfinebins);
% profiles_combine=combinelanes(profiles_align,numfinebins)
% Need to combine subdivisions for best estimate of scan down lanes.
%  Rhiju, August 31, 2003.
numpixels =size(profiles_align,1);
numlanes  =size(profiles_align,2);

%Now combine profiles within one lane together...
count=0;
%numfinebins
%numlanes
for k=1:numfinebins:numlanes
    count=count+1;
    profiles_combine(:,count) = sum( profiles_align( :, k + ...
						     [0:numfinebins-1]) ,2);
end

%figure(1)
%hold off;image(profiles_combine/(numfinebins*10));
%title(['Profiles down each lane, rescaled'])
