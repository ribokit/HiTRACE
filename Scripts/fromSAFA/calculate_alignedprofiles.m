function profiles_align = calculate_alignedprofiles(profiles,alignmentpoints,reflane);
numfinelanes = size(profiles,2);
numpixels = size(profiles,1);

numanchorlines = size(alignmentpoints,1);
%Pad beginning and end of profiles with zeros...
offset = numpixels;
profiles_offset = zeros(numpixels+2*offset, numfinelanes);
profiles_offset(offset+1:offset+numpixels,:) = profiles;

profiles_align = 0*profiles;

if numanchorlines ==0 profiles_align= profiles;return;end;

for k=1:numfinelanes
    % Just do a shift of the first block, above the first anchorline.
    beginpixel_new = 1; endpixel_new = floor(alignmentpoints(1,reflane)); 
    endpixel_old   = floor(alignmentpoints(1,k)); 
    beginpixel_old = endpixel_old + beginpixel_new - endpixel_new;
    profiles_align([beginpixel_new:endpixel_new],k) = ...
        profiles_offset([beginpixel_old:endpixel_old] + offset ,k);
    
    %Do rest of blocks, stretching to match sections between anchorlines.
    for j=1:numanchorlines-1
        beginpixel_new = floor(alignmentpoints(j,reflane))+1 ;
        endpixel_new   = floor(alignmentpoints(j+1,reflane));
        beginpixel_old = floor(alignmentpoints(j,k))+1 ;
        endpixel_old   = floor(alignmentpoints(j+1,k)) ;
        profile_temp =...
            profiles_offset([beginpixel_old:endpixel_old] + offset ,k);
        length_old = endpixel_old - beginpixel_old + 1;
        length_new = endpixel_new - beginpixel_new + 1;
        pixels_temp = (endpixel_new-beginpixel_new)*(0:length_old-1)/ (length_old-1) + beginpixel_new;
        %Use multiplicative factor, to preserve total counts on
        %gel (a "Jacobian" factor)?
        profile_stretch = (length_old/length_new)*interp1(pixels_temp, profile_temp,beginpixel_new:endpixel_new);
%        profile_stretch = interp1(pixels_temp, profile_temp,beginpixel_new:endpixel_new);
        profiles_align([beginpixel_new:endpixel_new],k) = profile_stretch';
    end 
end
