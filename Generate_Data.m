function [ Data, Covariates ] = Generate_Data()

beta_1 = 1; % starting point
beta_2 = 2; % slope

NumOfDataPoints = 20;

sigma = 0.5;

Covariates = [1:NumOfDataPoints]'; %vector ( 1 2 3 ... n )^T

Means = beta_1 + Covariates.*beta_2;

plot(Covariates, Means) % a line

Data = Means + sigma*randn(NumOfDataPoints,1);

hold on

scatter(Covariates, Data) % dots
end

