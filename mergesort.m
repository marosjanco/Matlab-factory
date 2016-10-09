function a = mergesort(a)
len = length(a);            % + c0

% T(1) = c

if len>1
   split = ceil(len/2);     % + c1
   
   a1 = a(1:split);
   a2 = a(split+1:end);     % + n*c2
   
   a1 = mergesort(a1);      
   a2 = mergesort(a2);      % + 2*T(n/2)
   
   i=1; j=1; pos=1;          % + c3
   
   % Merge STEP              % + c4*n + c5
       while i<=length(a1) && j<=length(a2)
           if a1(i)<a2(j)
               a(pos)=a1(i);
               i=i+1;
           else
               a(pos)=a2(j);
               j=j+1;
           end
           pos = pos + 1;
       end

       if i>length(a1)
           for k=j:length(a2)
               a(pos)=a2(k);
               pos = pos + 1;
           end
       else
           for k=i:length(a1)
               a(pos)=a1(k);
               pos = pos + 1;
           end  
       end
end
end

% z = [1,2,-4,0,0,15,188];
% mergesort(z)
    