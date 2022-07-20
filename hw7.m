clc
clear all
close all
%% 21.13
d13 = [0 2 4 6 8 10 12 14 16;0 0.7 1.8 3.4 5.1 6.3 7.3 8.0 8.4];
disp('문제 21.13')

% 중심차분
i = 6;
h = d13(1,i+1)-d13(1,i-1);
disp('중심차분')
fc = (d13(2,i+1)-d13(2,i-1))/(2*h)

% 전향차분
i = 6;
h = d13(1,i+1)-d13(1,i);
disp('전향차분')
ff = (d13(2,i+1)-d13(2,i))/h
% 후향차분
i = 6;
h = d13(1,i)-d13(1,i-1);
disp('후양차분')
fb = (d13(2,i)-d13(2,i-1))/h
%% 21.17
x17 = [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
fx17 = [0.05399 0.12952 0.24197 0.35207 0.39894 0.35207 0.24197 0.12952 0.05399];

% plot(x17,fx17)

disp('문제 21.17')
% f' 중심차분
for i = 2 : length(x17)-1
    h = x17(i+1)-x17(i-1);
    f17c(i-1) = (fx17(i+1)-fx17(i-1))/h;
end

% f''
for i = 2 : length(x17)-1
    h = (x17(i+1)-x17(i-1))/2;
    f17cd(i-1) = (fx17(i+1)-2*fx17(i)+fx17(i-1))/h^2;
end

% figure;
% plot(x17(2:length(x17)-1),f17c)
% figure;
% plot(x17(2:length(x17)-1),f17cd)

disp('변곡점은 f의 도함수 가 0 인 5번 데이터이다.')
%% 21.19
syms x
f19 =exp(-2*x)-x;
diff(f19,1);
disp('문제 21.19')
% x = 2
x=2;
disp('x=2 일때 도함수')
df = - 2*exp(-2*x) - 1

f = @(x) exp(-2*x)-x;
disp('dx 0.5 -> 0.001')

j = 1;
for h = 0.5:-0.001:0.01
    f19a(1,j) = (f(2+h)-f(2-h))/(2*h); % 중심
    f19a(2,j) = (f(2+h)-f(2))/h; % 전향
    f19a(3,j) = (f(2)-f(2-h))/h; % 후향
    j = j+1;
end
n = [0.5:-0.001:0.01];
figure;
plot(n,f19a(1,:))
hold on
plot(n,f19a(2,:))
plot(n,f19a(3,:))
plot(0,df, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
legend({'중심','전향','후향'})
hold off
%% 21.36
disp('문제 21.36')
L = 600;
E = 50000;
I = 30000;
w0 = 2.5;

th = @(x) w0/(120 * E * I * L) * (-5*x^4 + 6*L^2*x^2 - L^4);

% 수치적분을 이용하여 보의 처짐 계산
disp('보의 처짐')

h = 50;

LL = [0:h:L];

% simpson 1/3
simp3 = 0;
j = 2;
for i = 2:2:length(LL)
    if rem(L,2) == 0
        if i >length(LL)-1
            simp3(j) = simp3(j-1) + h * (th(LL(i-1))+th(LL(i)))/2;
        else
            simp3(j) = simp3(j-1) + h/3 * (th(LL(i-1))+4*th(LL(i))+th(LL(i+1)));
        end
    else
        simp3(j) = simp3(j-1) + h/3 * (th(LL(i-1))+4*th(LL(i))+th(LL(i+1)));
    end
    j = j+1;
end

simp3*10^-2
disp('m이다.')

%수치적분 모멘트와 전단계산

L = 6;
E = 200*10^6;
w0 = 2.5*10^3;
I = 0.0003;
th2 = @(x) w0/(120 * E * I * L) * (-5*x^4 + 6*L^2*x^2 - L^4);


h = 0.125;
d = [0:h:L];
for i = 1:length(d)
    fx(i) = th2(d(i));
end

% 중심차분 모멘트

for i = 1:length(fx)
    if i == 1,M(i) = (fx(i+1)-fx(i))/h *E*I;
    elseif i == length(fx), M(i) = (fx(i)-fx(i-1))/h *E*I;
    else M(i) = (fx(i+1)-fx(i-1))/(2*h) *E*I;
    end
end

% 중심차분 전단력
for i = 1:length(fx)
    if i == 1,V(i) = (M(i+1)-M(i))/h;
    elseif i == length(fx), V(i) = (M(i)-M(i-1))/h;
    else V(i) = (M(i+1)-M(i-1))/(2*h);
    end
end

figure;
hold on
plot(d,M)
plot(d,V)
hold off
legend({'th','moment','V'})
%% 21.41

v = @(t) 2*t /sqrt(1+t^2);

h1 = 0.5;
h2 = 0.25;
t=5;

d1 = (v(t+h1)-v(t-h1))/(2*h1);
d2 = (v(t+h2)-v(t-h2))/(2*h2);

disp('richardson')
D = 4/3 * d2 - 1/3 * d1

% 오차
disp('오차')
e1 = abs((D-d1)/D)
e2 = abs((D-d2)/D)