alpha = 0.2135;
beta  = 0.01709;
R  = 0.0820578;

L = 1:0.5:30;

A = zeros(length(L),6);
t = 0;
row = 0;

for T=40:5:50

    fprintf('T\n\n');
    
    t =t+1;
    row = 0;
    
    for P = L
row = row + 1;
row
a = (roots([1,-beta-(R*T)/P,alpha/P,-alpha*beta/P]))'

        
        A(row,1)= P;
        b = imag(a)== 0;
        c = a(b);
        if T == 40
            if row <=36 && row >= 17
                    A(row,2:4)=a;
            else
                        A(row,2)=real(c(1));

            end
            

        end 
        if T ==45
           A(row,5)=real(c(1));
        end
        if T == 50
           A(row,6)=real(c(1));
        
        end
    end
end


hold on;
t = title('V-P curves with different values of temperatures');
plot(A(:,2),A(:,1), '*-');
plot(A(17:36,2),A(17:36,1), '*-');
plot([A(17:36,3);A(17:36,4)],[A(17:36,1);A(17:36,1)], '*-');

plot(A(:,5),A(:,1), 'o-');
plot(A(:,6),A(:,1), '.-');
x = xlabel('V');
y = ylabel('P');
h = legend('T = 40 °K(first roots)','T = 40 °K(first roots as a part of 3 real root sln)','T = 40 °K(second and third roots)','T = 45 °K','T = 50 °K');
set(h,'FontSize',16);
set(t,'FontSize',15);
set(x,'FontSize',15);
set(y,'FontSize',15);
hold off
