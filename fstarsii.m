clf;
syms x m
hold on;
for p=1:9
    u = inline((4/pi)*symsum((1/((2*m-1))*sin((2*m-1)*pi*x/2)),m,1,p));
ezplot(u,[-2 2]);
end
f = @(x) (x<0).*(-1) + (x>0).*1;
ezplot(f,[-2,2]);
L = legend('p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','f(x)');
set(L, 'Position', [.76, .19, .1, .25]);
ylabel('$f(x)$  or $f_{2p-1}^{\star}(x)$ at particular $p$','Interpreter','LaTex')
xlabel('x');
title('$$f(x)$$ and $$f_{2p-1}^{\star}(x), p = 1 \rightarrow 9$$','interpreter','latex','fontsize',10)

