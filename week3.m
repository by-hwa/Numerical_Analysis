clc
clear all
close all

A = [3 -0.1 -0.2;1 7 -.3;3 -0.2 10]
B = [7.85;-19.3; 71.4]

[L, U] = lu(A)

[m,n] = size(A);
C = [A B]

for i = 1:m
    for j = (i+1):m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i))*C(i,:)
    end
end

X = zeros(1,n)';


for i = m:-1:1
    X(i,1) = (C(i,n+1) - C(i,1:n)*X)/C(i,i);
end

