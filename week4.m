clc
clear all
close all

%%
figure
ezplot('1*x^3 - 20*x^2 - 100*x + 2000',[-100, 100])
grid on

% matlab에서만 사용가능한 특이함수
Fx = @(x) 1*x^3 - 20*x^2 - 100*x + 2000;
%figure
%ezplot(Fx,[-100, 100])
Fx(1)
fzero(Fx,0) %구간 설정 가능

% 위는 다항식에서만 성립한다.
A = [1 -20 -100 2000];
polyval(A,1) %계수함수가 있을 때 대입한 결과를 알려줌
R = roots(A) % 함수가 갖는 근을 전부 알려줌

A1 = poly(R) % R근을 갖는 다항식의 계수를 보여줌 A를 보여준다.

%% LS(Line Searching): Inc(증분법)

xL = -50;
xU = 50;
N = 25;

dx = (xU-xL)/N;
x = [xL:dx:xU];
L = length(x);

k = 1;
j = 1;
for i = 1:L-1
    F(i) = A(1)*x(i)^3 + A(2)*x(i)^2 + A(3)*x(i) + A(4);
    F(i+1) = A(1)*x(i+1)^3 + A(2)*x(i+1)^2 + A(3)*x(i+1) + A(4);
    S = F(i) * F(i+1);
    if S < 0
        xr(k) = abs((x(i)+x(i+1))/2);
        k = k + 1;
    elseif S==0
        if F(i) == 0
            xt(j) = x(i);
            j = j+1;
        else
            xt(j) = x(i+1);
            j = j+1;
        end
    end
end

%% LineSearch : Bisection(2분법)

xL = -50;
xU = 50;
Es = 0.001;
N = log2((xU-xL)/Es);
N = ceil(N);

for i = 1:N
    xr(i) = (xU+xL)/2;
    FL = A(1)*xL^3 + A(2)*xL^2 + A(3)*xL + A(4);
    FR = A(1)*xr(i)^3 + A(2)*xr(i)^2 + A(3)*xr(i) + A(4);
    FU = A(1)*xU^3 + A(2)*xU^2 + A(3)*xU + A(4);
    SL = FL * FR;
    SU = FU * FR;
    if SL <0
        xU = xr(i);
    else
        xL = xr(i);
    end
end
