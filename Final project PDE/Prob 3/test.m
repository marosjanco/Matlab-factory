
n = 1;
z = 3;
while (n < z)
    n = n + 1;
    n
end

display(n);

n = 1;
while 1
    n = n + 1;
    n
    if (n == z)
        break;
    end
end

display(n);
