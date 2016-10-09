N = 13;
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,2*N);
for i = 1:N
   a(i) = 2^i;
end

for i = 1:N
   b(i) = 3*2^(i-1);
end

for i = 1:N
   c(2*i-1) = a(i);
   c(2*i) = b(i);
end
