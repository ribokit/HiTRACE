function [ area_peak_corrected ] = unsaturate( saturated_array, diluted_array, sd_cutoff )

%Corrects area_peak arrays for saturating bands.  
%Saturated array is the set of arrays measured at a full concentration; 
%diluted array is the same samples, in the same order, run at diluted version of the final sample.

%The diluted samples allow for saturating band quantitation.  REQUIRES THAT
%THE TWO SETS OF ARRAYS BE IN THE SAME ORDER (e.g. column 1 of
%saturated_array corresponds to the same sample as is in column 1 of
%diluted_array.  ALSO REQUIRES THAT AREA_PEAK ARRAYS BE LISTED WITH n ROWS IN 1
%COLUMN.

% sd_cutoff adjusts the number of std. devs a residual has to be from the
% mean to be rejected as a saturated band.  A higher sd_cutoff requires
% more deviation from the mean residual for a band to be excluded as an
% outlier.

if length(saturated_array) ~= length(diluted_array)
    error('different arrays','The saturated arrays do not equal the number of diluted arrays!')
end;

[num_rows, num_cols] = size(saturated_array);

for i = 1:num_cols;
ap_residual{i} = saturated_array(:,i) - diluted_array(:,i)  %holds residuals of area_peak for saturated - diluted samples
end;

for i = 1:num_cols;
sd_ap_res{i} = std(ap_residual{i}) %sd_ap_res holds standard deviation of the set of residuals for a pair of original/diluted samples
end;

devs = []
for i = 1:num_cols %goes through different folding conditions
    for j = 1:num_rows; %goes through the all of the peaks.
        devs(j,i) = ap_residual{i}(j) / sd_ap_res{i}; %devs holds the number of standard deviations above the mean a residual is.
    end;
end;


%a peak is not considered in the minimization algorithm if its residual
%between the two columns > 1.5 standard deviations
%a new version of area_peak will be assembled for the original/diluted
%samples, not considering peaks whose residuals are abnormal.

small_resi_sat = {}; %matrix holding area_peak values, in same order, for all points whose residuals are < 1.5 std. devs away from mean resid
small_resi_dil = {};
resi_counter = 1; %used to avoid empty spaces in small_resi matrix


%assembles new arrays in the small_resi variables; only containing peaks
%fewer than sd_cutoff std devs away from mean residual
for i = 1:num_cols;
    for j = 1:num_rows;
        if devs(j,i) <= sd_cutoff && devs(j,i) >= -1*(sd_cutoff)
            small_resi_sat{i}(resi_counter) = saturated_array(j,i);  
            small_resi_dil{i}(resi_counter) = diluted_array(j,i);
            resi_counter = resi_counter + 1;
        end;
    end;
    resi_counter = 1;
end;


% Solving for ideal scaling factor c via partial differentiation.


c_numerator = zeros(1,num_cols); %initializing
c_denominator = zeros(1,num_cols);
c_factor = {}; %will contain the ideal scaling factor

for i = 1:num_cols;
    for j = 1:length(small_resi_sat{i});
    c_numerator(1,i) = c_numerator(1,i) + small_resi_sat{i}(j)*small_resi_dil{i}(j);
    c_denominator(1,i) = c_denominator(1,i) + small_resi_dil{i}(j)*small_resi_dil{i}(j);
    end;
end;

for i = 1:num_cols;
    c_factor{i} = c_numerator(1,i) / c_denominator(1,i);
end;


area_peak_corrected = [];

%area_peak_corrected gets the original value for area_peak if it was deemed
%to be non-saturating; if the peak was deemed saturating, it gets the value
%measured in diluted array scaled by the constant c (calculated to make
%other peaks overlay with saturated array)
for i = 1:num_cols;
    for j = 1:num_rows;
        if devs(j,i) <= sd_cutoff && devs(j,i) >= -1*(sd_cutoff)
            area_peak_corrected(j,i) = saturated_array(j,i);
        elseif (devs(j,i) >= sd_cutoff) || (devs(j,i) <= -1*(sd_cutoff))
            area_peak_corrected(j,i) = c_factor{i}*diluted_array(j,i);
        end;
    end;
end;



end

