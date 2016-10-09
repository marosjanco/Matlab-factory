function [ MLE ] = Bernoulli_MLE( Prob, Num_Of_Samples )

% Run this function with increasing number of samples in each estimate
% As Num_Of_Samples get larger, the estimator converges in distribution to a normal

% For large Num_Of_Samples, the variance of the MLE -> 0


Num_Of_Estimates = 1000; % fixed at a large number to get good plots

MLE = zeros(1, Num_Of_Estimates);
T   = zeros(1, Num_Of_Estimates);
%WE MAKE ONE ESTIMATE FROM A GROUP OF SAMPLES (that create ONE DATA)

for i = 1:Num_Of_Estimates
    % Sample some data from Bernoulli distribution
    Data = binornd(1, Prob, 1, Num_Of_Samples); %NOTE: Data row vector will be always different, i.e. is random.

    % Calculating MLE
    MLE(i) = mean(Data);
    
    % Calculate alternative random variable
    T(i)   = sqrt(Num_Of_Samples)*(MLE(i) - Prob);
end


% Plot both the MLE and alternative rv (scaled MLE)
figure(1)
subplot(1,2,1)
% Plot MLE : converges to true value with variance -> 0
hist(MLE)
subplot(1,2,2)
% Plot scaled MLE : converges in distribution to N(0,1)
hist(T)

% Try plotting the empirical CDFs of these estimators

end

