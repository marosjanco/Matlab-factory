format long;
out = dlmread('out.txt','\t',1,1);

figure()
title(['The roots for the 81 quadratic equations given in the Mastery Section', char(10),'\Re(r_1)\geq \Re(r_2)'])
title(['The roots for the 81 quadratic equations given in the Mastery Section,', char(10), 'with $$\Re(r_1)\geq \Re(r_2)$$'],'interpreter','latex') 
title(['The roots for the 81 quadratic equations given in the Mastery Section,', char(10), 'where ' '$$r_1$$' ' and ' '$$r_2$$' ' were found by ''+sign'' and ''-sign'' quadratic formula respectively'],'interpreter','latex') 
hold on;
plot(out(:,1),out(:,2),'.b');
plot(out(:,3),out(:,4),'.g');
hold off;
legend('r1','r2');
xlabel('Re(w)')
ylabel('Im(w)')
xlim([0 1]); ylim([0 1]);
axis equal
