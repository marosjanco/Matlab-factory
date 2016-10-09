
function B = transfer(a,b,c)

    B = zeros(a*b*c,4);
    for n = 1 : a*b*c
        
        r1 = mod(n-1,a*b)+1;
        k = (n-r1)/(a*b) + 1;
        
        r2 = mod(r1-1,a)+1;
        j = (r1-r2)/a + 1;
        
        i = mod(r2-1,a) + 1;
        
        B(n,:) = [n, i,j,k];
    end

end
% 
% function B = transfer(a,b)
% 
%     B = zeros(a*b,3);
%     for pos = 1 : a*b
% 
%         r1 = mod(pos-1,a)+1;
%         j = (pos-r1)/a + 1;
%         
%         i = mod(r1-1,a) + 1;
%         
%         B(pos,:) = [pos, i,j];
%     end
