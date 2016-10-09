format long;
Val = zeros(1,10);
for p = 1:10
syms x m
fstar = ((4/pi)*symsum((1/((2*m-1))*sin((2*m-1)*pi*(1/x)/2)),m,1,x));  
Val(p) = eval(subs(fstar, x, p*100))-1; %substituing function fstar at x = (10*p) and then evaluating the symbolic expression
end 
disp(Val);