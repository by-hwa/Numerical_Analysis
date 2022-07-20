clc
clear all
close all
%% 함수
syms x
f = 20 * (x-3)*(x-7)*(x+1)

expand(f)
ezplot(f)

%% 적분 구간
a = -10;
b = 10

%% 오직 매트랩

F = int(f)
Fa = subs(F,x,a);
Fb = subs(F,x,b);

Int_Mat = Fb - Fa
int_mat = int(f,[a,b])

%% 데이터
dx = 1;
D = [a:dx:b];
L = length(D);
N = L-1;

for i = 1:L
    fD(i) = 20*D(i)^3 - 180*D(i)^2 + 220*D(i) + 420;
end

figure
plot(D,fD, 'ko')

%% 합성 사다리꼴
h = (b-a)/N;
Int_trap = 0;
for i = 1:L
    if i == 1
        Int_trap(i) = (h/2) * fD(i);
    elseif i == L
        Int_trap(i) = Int_trap(i-1) + (h/2) * fD(i);
    else
        Int_trap(i) = Int_trap(i-1) + (h/2) * fD(i)*2;
    end
end
Int_trap(end)

%실시간 적분
% 두데이터의 평균
Int_trapR = 0;
for i =2:L
    Int_trapR(i) = Int_trapR(i-1) + h * (fD(i-1)+fD(i))/2;
end
Int_trapR(end)

% 직사각형
% 실시간 데이터 값의 크기

Int_recR = 0;
for i =2:L
    Int_recR(i) = Int_recR(i-1) + h * fD(i);
end
Int_recR(end)

%% 데이터 romberg

dx = (b-a)/8;
D = [a:dx:b];
L = length(D);
N = L-1;

fD = 0;

for i = 1:L
    fD(i) = 20*D(i)^3 - 180*D(i)^2 + 220*D(i) + 420;
end

figure
plot(D,fD, 'ko')

%% Romberg
h = (b-a)/1;
Int_trap1 = 0;
for i =1:L-N
    Int_trap1 = Int_trap1 + h * (fD(i+N)+fD(i))/2;
end
Int_trap1

h = (b-a)/2;
Int_trap2 = 0;
for i =1:N/2:L-N/2
    Int_trap2 = Int_trap2 + h * (fD(i+N/2)+fD(i))/2;
end
Int_trap2

h = (b-a)/4;
Int_trap3 = 0;
for i =1:N/4:L-N/4
    Int_trap3 = Int_trap3 + h * (fD(i+N/4)+fD(i))/2;
end
Int_trap3

h = (b-a)/8;
Int_trap4 = 0;
for i =1:N/8:L-N/8
    Int_trap4 = Int_trap4 + h * (fD(i+N/8)+fD(i))/2;
end
Int_trap4

Int_R = [Int_trap1;Int_trap2;Int_trap3;Int_trap4]
L_R = 4;
for k = 2:L_R
    for j = 1:L_R-k+1
        Int_R(j,k) = (4^(k-1)*Int_R(j+1,k-1) - Int_R(j,k-1))/(4^(k-1)-1);
    end
end
Int_R

%% romberg 실시간.

Int_R = 0;
k = 2;

for i = 2:L
    if rem(i,2) == 1
        Int_T1 = 2*h*(fD(i)+fD(i-2))/2;
        Int_T2 = h*(fD(i)+fD(i-1))/2 + h*(fD(i-1)+fD(i-2))/2;
        Int_R(k) = Int_R(k-1) + (4/3)*Int_T2 - (1/3)*Int_T1;
        k = k +1;
        i
    end
end
Int_R