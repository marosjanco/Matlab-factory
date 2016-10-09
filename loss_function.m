x = 0:0.01:1;
P = x.^4.*(1-x).^4;
n = 10;
close all
hold all
for ind=2:n-2
    argmin = ind/n;
    L = (x-argmin).^2;
    fprintf(num2str(argmin),'\n');
%     plot(x,L)
%     plot(x,P)
    plot(x,L.*P,'DisplayName',num2str(argmin))
%     legend(gca,'show')
%     legend('loss','probability','multiplied')
%     pause(2)
end
legend(gca,'show')
