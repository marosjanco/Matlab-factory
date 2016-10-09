a = [-5 9 -7 2];
a = [3 -39 -47 210];

a = flip(a);
b = [(4*a(1)*a(3) - a(2)^2-a(1)*a(4)^2 ) (a(2)*a(4)-4*a(1)) -a(3)];
b = flip(b);
r = roots([1,b]);
r = max(real(r));

c1 = [1, (a(4)/2 + sqrt(a(4)^2/4+r-a(3))), (r/2 + sqrt(r^2/4-a(1)))];
c2 = [1, (a(4)/2 - sqrt(a(4)^2/4+r-a(3))), (r/2 - sqrt(r^2/4-a(1)))];
roots(c1);
roots(c2);

c1 = [1, (a(4)/2 + sqrt(a(4)^2/4+r-a(3))), (r/2 - sqrt(r^2/4-a(1)))];
c2 = [1, (a(4)/2 - sqrt(a(4)^2/4+r-a(3))), (r/2 + sqrt(r^2/4-a(1)))];
roots(c1);
roots(c2);
R = [
    %-7,-3,2,5;
    1+1i, 1-1i, 2i -2i;
    7,7,7,7;
    1,1,1,2;
    3,3,3,1;
    -5,-5,1+1i,1-1i;
    -4, 3, 1+1i, 1-1i;
    10, 10, -110 -110;
    10, 10, 1, -3;
    -3, -3, 10, 1;
    2, 1, -1, -3
    ];
R = [
    4, 2+8*1i, 2-8*1i;
    -7, -7, -7;
    8, 8, 8;
    0, 0, 5;
    4, 4, 0;
    4, 4, 5;
    1, 2, 0;
    0, 2, 3;
    1, 0, 3;
    1, 2, 3;
    1, -18, 300
    ];
R = [
    16 + 1i, 16 - 1i;
    -16*1i,16*1i;
    1 + 16*1i, 1 - 16*1i;
    2, 2;
    -3, -3;
    16,32;
    16,-32;
    -16,32;
    -16,-32

    ];
S = size(R);
C = zeros(length(R),S(2)+1);
RR = zeros(S);
for i=1:length(R)
    C(i,:) = poly(R(i,:));
    RR(i,:) = (roots(C(i,:)))';

end