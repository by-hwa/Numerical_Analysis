clc
clear all
close all

%% Function
syms x
f = 20*(x-3)*(x-7)*(x+1);
f = expand(f)
ezplot(f)

%% matlab
f_dot = diff(f)
f_ddot = diff(diff(f))

a = 5;

Dif_Mat = subs(f_dot, x, a)
Dif_Mat = double(Dif_Mat)

%% Data
dx = 1;
D = [0:dx:8];
L = length(D);
N = L-1;
C = coeffs(f);
C = double(C);

for i = 1:L
    fD(i) = C(1) + C(2)*D(i) + C(3)*D(i)^2 + C(4)*D(i)^3;
end

figure;
plot(D,fD, 'ko')
hold on
ezplot(f, [0,8])
%% 수치 미분
% Forward
Dif_F1 = (fD(a/dx + 2) - fD(a/dx + 1))/dx
Dif_F2 = (-fD(a/dx + 3) + 4*fD(a/dx + 2) - 3*fD(a/dx + 1))/(2*dx)

% Backward
Dif_B1 = (fD(a/dx + 1) - fD(a/dx))/dx
Dif_B2 = (3*fD(a/dx + 1) - 4*fD(a/dx) + fD(a/dx - 1))/(2*dx)

% Centerd
Dif_C1 = (fD(a/dx + 2) - fD(a/dx))/(2*dx)
Dif_C2 = (-fD(a/dx + 3) + 8*fD(a/dx + 2) - 8*fD(a/dx) + fD(a/dx - 1))/(12*dx)

E_F1 = Dif_Mat - Dif_F1
E_F2 = Dif_Mat - Dif_F2

E_B1 = Dif_Mat - Dif_B1
E_B2 = Dif_Mat - Dif_B2

E_C1 = Dif_Mat - Dif_C1
E_C2 = Dif_Mat - Dif_C2

%% Real
Dif_B1 = 0;
Dif_B2 = 0;
Dif_B3 = 0;

for i = 2:L
    Dif_B1(i) = (fD(i)-fD(i-1))/dx;
end

for i = 3:L
    Dif_B2(i) = (3*fD(i) - 4*fD(i-1) + fD(i-2))/(2*dx)
end

figure;
ezplot(f_dot,[0,8])
hold on
plot(D,Dif_B1,'ko')
plot(D,Dif_B2, 'ro')

%% Romberg
Dif_B11(1) = 0;
Dif_B11(2) = 0;
Dif_B12(1) = 0;
Dif_B12(2) = 0;
Dif_R(1) = 0;
Dif_R(2) = 0;

for i = 3:L
    Dif_B11(i) = (fD(i)-fD(i-1))/dx;
    Dif_B12(i) = (fD(i)-fD(i-2))/(2*dx);
    Dif_R(i) = (4/3)*Dif_B11(i) - (1/3)*Dif_B12(i);
end

plot(D,Dif_R,'mo')

Dif_B11(1) = 0;
Dif_B11(2) = 0;
Dif_B11(3) = 0;
Dif_B11(4) = 0;
Dif_B12(1) = 0;
Dif_B12(2) = 0;
Dif_B12(3) = 0;
Dif_B12(4) = 0;
Dif_R(1) = 0;
Dif_R(2) = 0;
Dif_R(3) = 0;
Dif_R(4) = 0;

for i = 5:L
    Dif_B11(i) = (3*fD(i) - 4*fD(i-1) + fD(i-2))/(2*dx);
    Dif_B12(i) = (3*fD(i) - 4*fD(i-2) + fD(i-4))/(4*dx)
    Dif_R2(i) = (4/3)*Dif_B11(i) - (1/3)*Dif_B12(i);
end

plot(D,Dif_R2,'go')