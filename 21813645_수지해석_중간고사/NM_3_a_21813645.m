clc
clear all
close all

k = 10;
m1 = 2;
m2 = 3;
m3 = 2.5;
g = 9.81;

k1 = k;
k2 = k;
k3 = k;


K = [k1+k2 -k2 0;-k2 k2+k3 -k3;0 -k3 k3];
mg = [m1*g; m2*g; m3*g;];

disp('x1,x2,x3ÀÇ º¯À§')
X = inv(K) * mg