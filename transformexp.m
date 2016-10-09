


n = 5000;
mu = 1;
r1 = unifrnd(0,1,1,n); %standard uniform (0,1), 1xn row vector
r2 = -(1/mu)*log(ones(1,n)-r1);


figure(1)
subplot(1,2,1)
hold on;
scatter(0*ones(1,n),r1);
scatter(ones(1,n),r2);
hold off;
subplot(1,2,2)
hist(r2,20)
