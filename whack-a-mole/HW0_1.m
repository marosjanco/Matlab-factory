plot(1:5,1:5)
i = 1;
j = 2;
length_hit = 10;
k = 4;

text = {['Hitting the true positions in order']; ...
                ['Coordinates hit = (' num2str(i) ',' ...
                num2str(j) '); ' num2str(length_hit-k) ' left)']}
title(text)