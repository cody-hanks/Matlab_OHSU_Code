clc;
clear;
close all;

a = 10*randn(1000,1);
b = randn(1000,1);
plot(a,b,'.');
%R = randi(100,2,2)*.001;
R = [cosd(75) -sind(75); sind(75) cosd(75)];
G= R*[a,b]';
x=G(1,:);
y=G(2,:);
hold on; 
plot(x,y,'G.');
[coeff]=pca([x;y]');
coeff
R
G1=coeff*[x;y];
hold on;
plot(G1(1,:),G1(2,:),'r.')