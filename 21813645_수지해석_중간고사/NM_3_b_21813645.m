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

C = [K mg];



for i = 1:3
    for j = i+1:3
    C(j,:) = C(j,:)-(C(j,i)/C(i,i))*C(i,:);
    end
end

x = zeros(3,1);

for i = 3: -1:1
    x(i,:) = (C(i,4)-(C(i,1:3)*x))/C(i,i);
end


x

