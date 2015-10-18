clear
clc
mex  -g Harris.cpp
fprintf('ok\n');

I =imread('1.jpg');
if(size(I,3)==3)
    I=rgb2gray(I);
end
I=double(I);
[x y]=harris(I,2,20,12);
figure;
imshow(uint8(I));
hold on
plot(x,y,'+');

