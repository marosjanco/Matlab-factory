year = [2014	2013	2012];
tot_asset = [2482643	1062037	1316664];
liab = [1285087	819256	308581];

tot_asset = [2423119	1823016];
liab = [3441964	2038689];

tot_asset = [1426683	1257845	850628];
liab = [291707	261288	214409];

x = year(1:end);
y = tot_asset - liab

p = polyfit(x,y,1);
f = @(x) p(1)*x+p(2);



hold on;
plot(x,y,'o');
ezplot( f,[min(x) 2016] );
plot(2016,f(2016),'*');

title('Linear regression estimate for 2016 total equity');
xlabel('Year');
ylabel('Total equity');
xlim([min(x) 2016]);
ylim([min(y) max(y)]);
legend('known total equity','linear fit for the known data','approximated total equity for 2016');
grid on;
hold off;
display(f(2016));