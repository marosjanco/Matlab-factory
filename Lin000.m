format compact;

 m = 10;
 n = 20;

 
 a = randn(m,n);
 b = randn(m,n)<-0.5;
 
 sum(sum(a.*b))
 
 trace(a'*b)