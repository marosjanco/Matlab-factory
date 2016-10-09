clf; hold on;
for p=1:20
   fplot(s(2*p-1),[-2 2]);
end
f = @(x) (x<0).*(-1) + (x>0).*1;
ezplot(f,[-2,2]);
L = legend('p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','f(x)');
set(L, 'Position', [.76, .19, .1, .25]);
ylabel('$f(x)$  or $s_{2p-1}(x)$ at particular $p$','Interpreter','LaTex')
xlabel('x');
title('$$f(x)$$ and $$s_{2p-1}(x), p = 1 \rightarrow 9$$','interpreter','latex','fontsize',10)

