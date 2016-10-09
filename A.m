N = 32;
S = zeros(N-1,N-1);

for i = 1: N-1
    for j = 1: N-1
        S(i,j) = sin(i*j*pi/N);
    end
end

y = (1:N-1)';
true_res = S*y;
display(true_res);

SS = zeros(N/2-1,N/2-1);

for i = 1: N/2-1
    for j = 1: N/2-1
        SS(i,j) = sin(i*j*pi/(N/2));
    end
end

T = zeros(N/2,N/2);

for i = 1: N/2
    for j = 1: N/2
        T(i,j) = sin((2*i-1)*j*pi/N);
    end
end

% display(SS);
% display(T);

P = [[1,0,0,0,0,0,-1],
 [0,1,0,0,0,-1,0],
 [0,0,1,0,-1,0,0],
 [1,0,0,0,0,0,1],
 [0,1,0,0,0,1,0],
 [0,0,1,0,1,0,0],
 [0,0,0,1,0,0,0]];

a_s = P*y;  % a is less
a = a_s(1:N/2-1);
s = a_s(N/2:N-1);

display(a);
display(s);

display(SS*a);
display(T*s);




T = zeros(2,2);

for i = 1: 2
    for j = 1: 2
        T(i,j) = sin((2*i-1)*j*pi/2);
    end
end
