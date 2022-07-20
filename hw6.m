%% 19.3
clc
clear all
close all

disp('19.3')
a = 0;
b = pi/2;
dx = 0.1;
D = [a:dx:b];

for i = 1:length(D)
    FD(i) = 8 + 4*cos(D(i));
end

% 메트렙
syms x
f = 8 + 4*cos(x);
disp('해석적 방법')
int(f,[a,b])
sol = 4*pi + 4

% 합성사다리꼴
h = (b-a)/length(D);
trap = 0;

for i = 1:length(D)
    if i == 1
        trap(i) = h/2 * FD(i);
    elseif i == length(D)
        trap(i) = trap(i-1)+ h/2 * FD(i);
    else
        trap(i) = trap(i-1) + h/2 * FD(i)*2;
    end
end

disp('합성사다리꼴')
disp('n = 16')
trap(end)
disp('오차')
abs((sol-trap(end))/sol)
L = length(D);


% 합성 simpson 1/3
simp3 = 0;
j = 2;
4*pi + 4
for i = 2:2:L
    if rem(L,2) == 0
        if i == L
            simp3(j) = simp3(j-1) + h * (FD(i-1)+FD(i))/2;
        else
            simp3(j) = simp3(j-1) + h/3 * (FD(i-1)+4*FD(i)+FD(i+1));
        end
    else
        simp3(j) = simp3(j-1) + h/3 * (FD(i-1)+4*FD(i)+FD(i+1));
    end
    j = j+1;
end

disp('합성 sipson 1/3')
simp3(end)
disp('오차')
abs((sol-simp3(end))/sol)

% simpson 1/3
dx = b/3;
D = [a:dx:b];

for i = 1:length(D)
    fd(i) = 8 + 4*cos(D(i));
end


h = (b-a)/2;
disp('simpson 1/3')
I = 1/3*h*(fd(1)+4*fd(2)+fd(3))
disp('오차')
abs((sol-I)/sol)

% simpson 3/8
dx = b/4;
D = [a:dx:b];

for i = 1:length(D)
    fd(i) = 8 + 4*cos(D(i));
end


h = (b-a)/3;
disp('simpson 3/8')
I = 3/8*h*(fd(1)+3*fd(2)+3*fd(3)+fd(4))
disp('오차')
abs((sol-I)/sol)
%% 19.5
disp('19.5')
a = 0;
b = 1.2;
dx = 0.1;
D5 = [a:dx:b];

syms x
f5 = exp(-x);

int(f5,[a,b]);

for i = 1:length(D5)
    fd5(i) = exp(-D5(i));
end

disp('해석적 방법')
sol = 1 - exp(-6/5)

% 사다리꼴 공식
h = (b-a)/length(D5);

trap5 = 0;
for i = 2:length(D5)
    trap5(i) = trap5(i-1) + h * (fd5(i-1)+fd5(i))/2;
end

disp('합성 사다리꼴 방식')
trap5(end)
disp('오차')
abs((sol-trap5(end))/sol)

% simpson3
simp = 0;
j = 2;
for i =2:2:length(D5)
   simp(j) = simp(j-1) + h * 1/3 * (fd5(i-1)+4*fd5(i)+fd5(i+1));
   j = j+1;
end

disp('합성 simpson1/3 방식')
simp(end)
disp('오차')
abs((sol-simp(end))/sol)
%% 19.25
disp('19.25')
d25 = [0 50 100 150 225 300 375 450 600;0 30 40 40 50 50 60 80 100];

[m, n] = size(d25);


% 사다리꼴
trap25 = 0;
for i = 2:n
    trap25(i) = trap25(i-1) + (d25(1,i)-d25(1,i-1))*((d25(2,i)+d25(2,i-1))/2);
end
disp('합성 사다리꼴')
trap25(end)

% simpson 1/3
simp3 = 0;
j = 2;
for i = 2:2:n
    simp3(j) = simp3(j-1) + (d25(1,i+1)-d25(1,i-1))*((4*d25(2,i)+d25(2,i-1)+d25(2,i+1))/6);
    j = j+1;
end

disp('합성 simpson 1/3')
simp3(end)

% simpson 3/8
simp38 = 0;
j = 2;
for i = 2:3:n
    if n-i<4
        simp38(j) = simp38(j-1) + (d25(1,i+1)-d25(1,i-1))*((4*d25(2,i)+d25(2,i-1)+d25(2,i+1))/6);
    else
        simp38(j) = simp38(j-1) + (d25(1,i+2)-d25(1,i-1))*((3*d25(2,i)+d25(2,i-1)+3*d25(2,i+1)+d25(2,i+2))/8);
    end
    j = j+1;
end

disp('simpson 3/8 + 1/3')
simp38(end)
%% 19.26
disp('19.26')
d26 = [0 4 8 12 16 20 24 28 30;0 18 31 42 50 56 61 65 70];
[m, n] = size(d26);

% 사다리꼴
trap26 = 0;
for i = 2:n
    trap26(i) = trap26(i-1) + (d26(1,i)-d26(1,i-1))*((d26(2,i)+d26(2,i-1))/2);
end
disp('합성 사다리꼴')
trap26(end)
disp('평균속도')
trap26(end)/30
% simpson 1/3
simp3 = 0;
j = 2;
for i = 2:2:n
    simp3(j) = simp3(j-1) + (d26(1,i+1)-d26(1,i-1))*((4*d26(2,i)+d26(2,i-1)+d26(2,i+1))/6);
    j = j+1;
end

disp('합성 simpson 1/3')
simp3(end)
disp('평균속도')
simp3(end)/30
% simpson 3/8
simp38 = 0;
j = 2;
for i = 2:3:n
    if n-i<4
        simp38(j) = simp38(j-1) + (d26(1,i+1)-d26(1,i-1))*((4*d26(2,i)+d26(2,i-1)+d26(2,i+1))/6);
    else
        simp38(j) = simp38(j-1) + (d26(1,i+2)-d26(1,i-1))*((3*d26(2,i)+d26(2,i-1)+3*d26(2,i+1)+d26(2,i+2))/8);
    end
    j = j+1;
end
disp('simpson 3/8 + 1/3')
simp38(end)
disp('평균속도')
simp38(end)/30