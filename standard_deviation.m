function [ st_dev ] = standard_deviation( data )
format long   % expressing as a fixed-decimal format with 15 digits after the decimal point for double values
n = length(data); % number of observations
average = (1/n)*sum(data); % mean of the observations
st_dev = sqrt((1/(n-1))*sum((data - average*ones(1,n)).^2)); % emperical standard deviation
end

