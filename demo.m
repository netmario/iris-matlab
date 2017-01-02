close all;
[iris1,mask1] = iris('eye.jpg');
[iris2,mask2] = iris('eye2.jpg');

%save('codes.mat','iris1','iris2','mask1','mask2')
%load('codes.mat')

dist = smallest_distance(iris1,iris2,mask1,mask2)
