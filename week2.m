clc
clear all
close all

A = [1 4 2;5 7 3;1 3 9]
B = [3 4;2 4;1 6]
C = [2 3 22;2 45 32;231 54 2]

D = inv(A)

A
A' 

det(A)
det(C)

eig(A) % 고유치 모든 고유치의 곲은 det(a)이다


[EV, E] = eig(A) % 고유벡터 고유치 따로 구분가능


A = [1 4 2;5 7 3;1 3 9]
B = [3;4;1]

X = inv(A)*B