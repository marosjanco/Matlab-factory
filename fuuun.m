hold all;
for ind=1:3
x=[0:0.1:10];
plot(x, sin(x)+ind, 'DisplayName',['sin + ' num2str(ind)]);
end
legend(gca,'show')