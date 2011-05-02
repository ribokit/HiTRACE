function Z=hSHAPE(area_in)
%This calculates the SHAPE reactivity for the inputted data 

data_size=length(area_in);
data_number=size(area_in,2);

for i=1:data_number
    buffer=[];
    buffer=area_in(:,i);
Q=quartiles(buffer);

IQR=Q(3)-Q(1);
max_Q=IQR*1.5+Q(3);

% determine extreme Q3 outliers (e.g., x > Q1 + 3*IQR)
iz = find(buffer>Q(1)+1.5*IQR);

%Remove outliers
Z_post=buffer;
Z_post(iz)=Q(2);
Z_final=sort(Z_post);

%identify the high 10% and low 10%
per=.9;
grouping=floor(data_size*(1-per));

hrange(1)=mean(Z_final(end-grouping:end));
hrange(2)=mean(Z_final(1:grouping));

%Z=(area_in-hrange(2))./(hrange(1)-hrange(2));
Z(:,i)=buffer./(hrange(1));

end

function Q=quartiles(X)

Nx = size(X,1);

% compute mean
mx = mean(X);

% compute the standard deviation
sigma = std(X);

% compute the median
medianx = median(X);

% STEP 1 - rank the data
y = sort(X);

% compute 25th percentile (first quartile)
Q(1) = median(y(find(y<median(y))));

% compute 50th percentile (second quartile)
Q(2) = median(y);

% compute 75th percentile (third quartile)
Q(3) = median(y(find(y>median(y))));