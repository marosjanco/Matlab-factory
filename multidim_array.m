mat1 = [1/sqrt(2),-sqrt(3)/2;sqrt(3)/2,1/sqrt(2)];
mat2 = [-1,0;0,1];
mat3 = [1,1/2;0,1];
mat4 = [2,0;0,1/2];

matArray(:,:,1) = mat1;
matArray(:,:,2) = mat2;
matArray(:,:,3) = mat3;
matArray(:,:,4) = mat4;

fPoints = [0,0;0,2;1,2;0,2;0,1;0.75,1]'

matNumber = input('1: rotation; 2: reflection; 3: shear; 4: stretch-shrink. Choose? ');

plot(fPoints(1,:),fPoints(2,:),'r');
axis equal;
axis([-2 2 -1 3]);

% if matNumber == 1
% mat = mat1;
% elseif matNumber == 2
% mat = mat2;
% elseif matNumber == 3
% mat = mat3;
% else
% mat = mat4;
% end
% 


% matArray = {mat1, mat2, mat3, mat4};

% transformedFPoints = matArray{matNumber}*fPoints;

transformedFPoints = matArray(:,:,matNumber)*fPoints;
hold on;
plot(transformedFPoints(1,:),transformedFPoints(2,:),'b');
hold off;