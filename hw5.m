clc
clear all
close all
%% 9-6 해 탐색 법
x1 = [0];
x2 = [0];
x3 = [0];
i = 2;
% gauss seidal
while 1
    x1(i) = (27-(2*x2(end))-x3(end))/10;
    x2(i) = (-61.5+(3*x1(end))-(2*x3(end)))/-3;
    x3(i) = (-21.5-x1(end)-x2(end))/6;
    
    E1 = abs((x1(i)-x1(i-1))/x1(i))*100;
    E2 = abs((x2(i)-x2(i-1))/x2(i))*100;
    E3 = abs((x3(i)-x3(i-1))/x3(i))*100;
    
    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002, break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 gauss seidal 법으로 찾은 x1, x2 x3의 값이다')
x1(end)
x2(end)
x3(end)
disp('시행횟수 '), i


x1 = [0];
x2 = [0];
x3 = [0];
i = 2;

% jacobi
while 1
    x1(i) = (27-(2*x2(i-1))-x3(i-1))/10;
    x2(i) = (-61.5+(3*x1(i-1))-(2*x3(i-1)))/-3;
    x3(i) = (-21.5-x1(i-1)-x2(i-1))/6;
    
    E1 = abs((x1(i)-x1(i-1))/x1(i))*100;
    E2 = abs((x2(i)-x2(i-1))/x2(i))*100;
    E3 = abs((x3(i)-x3(i-1))/x3(i))*100;

    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002, break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 jacobi 법으로 찾은 x1, x2 x3의 값이다')
x1(end)
x2(end)
x3(end)
disp('시행횟수 '), i

%% 9-15 해 탐색 법
Ax = [0];, Ay = [0];, Ey = [0];, AB = [0];, BC = [0];
AD = [0];, BD = [0];, CD = [0];, DE = [0];, CE = [0];
i =  2;
ramb = 0.5;


while 1
    AB(i) = -0.8*BD(end);,AB(i) = (AB(i) * ramb) + (1-ramb) * AB(i-1);
    BC(i) = -74-(0.6*BD(end));,BC(i) = (BC(i) * ramb) + (1-ramb) * BC(i-1);
    AD(i) = DE(end) - (0.6*BD(end));,AD(i) = (AD(i) * ramb) + (1-ramb) * AD(i-1);
    BD(i) = -CD(end)/0.8;,BD(i) = (BD(i) * ramb) + (1-ramb) * BD(i-1);
    CD(i) = -24-(0.8*CE(end));,CD(i) = (CD(i) * ramb) + (1-ramb) * CD(i-1);
    DE(i) = -0.6*CE(end);,DE(i) = (DE(i) * ramb) + (1-ramb) * DE(i-1);
    CE(i) = BC(end)/0.6;,CE(i) = (CE(i) * ramb) + (1-ramb) * CE(i-1);
    Ax(i) = -AD(end);,Ax(i) = (Ax(i) * ramb) + (1-ramb) * Ax(i-1);
    Ay(i) = -AB(end);,Ay(i) = (Ay(i) * ramb) + (1-ramb) * Ay(i-1);
    Ey(i) = -0.8*CE(end);,Ey(i) = (Ey(i) * ramb) + (1-ramb) * Ey(i-1);

    E1 = abs((Ax(i)-Ax(i-1))/Ax(i))*100;
    E2 = abs((Ay(i)-Ay(i-1))/Ay(i))*100;
    E3 = abs((Ey(i)-Ey(i-1))/Ey(i))*100;
    E4 = abs((AB(i)-AB(i-1))/AB(i))*100;
    E5 = abs((BC(i)-BC(i-1))/BC(i))*100;
    E6 = abs((AD(i)-AD(i-1))/AD(i))*100;
    E7 = abs((BD(i)-BD(i-1))/BD(i))*100;
    E8 = abs((CD(i)-CD(i-1))/CD(i))*100;
    E9 = abs((DE(i)-DE(i-1))/DE(i))*100;
    E10 = abs((CE(i)-CE(i-1))/CE(i))*100;
    
    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002 && E4 < 0.002 && E5 < 0.002 && E6 < 0.002 && E7 < 0.002 && E8 < 0.002 && E9 < 0.002 && E10 < 0.002 , break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 Gauss-seidal 법으로 찾은 값이다')
    AB(end), BC(end), AD(end), BD(end), CD(end), DE(end), CE(end), Ax(end), Ay(end),Ey(end)
disp('시행횟수 '), i

Ax = [0];, Ay = [0];, Ey = [0];, AB = [0];, BC = [0];
AD = [0];, BD = [0];, CD = [0];, DE = [0];, CE = [0];
i =  2;
ramb = 0.5;


while 1
    AB(i) = -0.8*BD(i-1);
    BC(i) = -74-(0.6*BD(i-1));
    AD(i) = DE(i-1) - (0.6*BD(i-1));
    BD(i) = -CD(i-1)/0.8;
    CD(i) = -24-(0.8*CE(i-1));
    DE(i) = -0.6*CE(i-1);
    CE(i) = BC(i-1)/0.6;
    Ax(i) = -AD(i-1);
    Ay(i) = -AB(i-1);
    Ey(i) = -0.8*CE(i-1);

    AB(i) = (AB(i) * ramb) + (1-ramb) * AB(i-1);
    BC(i) = (BC(i) * ramb) + (1-ramb) * BC(i-1);
    AD(i) = (AD(i) * ramb) + (1-ramb) * AD(i-1);
    BD(i) = (BD(i) * ramb) + (1-ramb) * BD(i-1);
    CD(i) = (CD(i) * ramb) + (1-ramb) * CD(i-1);
    DE(i) = (DE(i) * ramb) + (1-ramb) * DE(i-1);
    CE(i) = (CE(i) * ramb) + (1-ramb) * CE(i-1);
    Ax(i) = (Ax(i) * ramb) + (1-ramb) * Ax(i-1);
    Ay(i) = (Ay(i) * ramb) + (1-ramb) * Ay(i-1);
    Ey(i) = (Ey(i) * ramb) + (1-ramb) * Ey(i-1);
    
    E1 = abs((Ax(i)-Ax(i-1))/Ax(i))*100;
    E2 = abs((Ay(i)-Ay(i-1))/Ay(i))*100;
    E3 = abs((Ey(i)-Ey(i-1))/Ey(i))*100;
    E4 = abs((AB(i)-AB(i-1))/AB(i))*100;
    E5 = abs((BC(i)-BC(i-1))/BC(i))*100;
    E6 = abs((AD(i)-AD(i-1))/AD(i))*100;
    E7 = abs((BD(i)-BD(i-1))/BD(i))*100;
    E8 = abs((CD(i)-CD(i-1))/CD(i))*100;
    E9 = abs((DE(i)-DE(i-1))/DE(i))*100;
    E10 = abs((CE(i)-CE(i-1))/CE(i))*100;
    
    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002 && E4 < 0.002 && E5 < 0.002 && E6 < 0.002 && E7 < 0.002 && E8 < 0.002 && E9 < 0.002 && E10 < 0.002 , break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 jacobi 법으로 찾은 값이다')
    AB(end), BC(end), AD(end), BD(end), CD(end), DE(end), CE(end), Ax(end), Ay(end),Ey(end)
disp('시행횟수 '), i

disp('연속대입법만으로 값을 찾을때, 발산하여 이완법을 적용하였습니다.')

%% 12-6 선형 연립 방정식


A = [2 -6 -1;-3 -1 7;-8 1 -2];
B = [-38;-34;-20];
% 역행렬
disp('역행렬')
X = inv(A) * B

% 소거법
C = [A B];
[m,n] = size(A);

for i = 1:m
    for j = (i+1) :m
        C(j,:) = C(j,:) - ((C(j,i)/C(i,i))*C(i,:));
    end
end

X = zeros(n,1);

for i = n:-1:1
    X(i,:) = (C(i,end) - (C(i,1:n)*X))/C(i,i);
end

disp('gauss 소거/대입법')
X

A = [A(3,:);A(1,:);A(2,:)];
B = [B(3,:);B(1,:);B(2,:)];

[L,U] = lu(A);

disp('LU분해법 Ax=B, LUx = B, LY=B, UX=Y')

Y = zeros(n,1);

for i = 1:m
    Y(i,:) = (B(i,:)-(L(i,1:3)*Y))/L(i,i);
end


X = zeros(n,1);

for i = m:-1:1
    X(i,:) = (Y(i,:) - (U(i,1:3)*X))/U(i,i);
end
X

% 탐색법
x1 = [0];
x2 = [0];
x3 = [0];

i = 2;
% gauss seidal
ramb = 1.2;

while 1
    x1(i) = (-20-x2(end)+2*x3(end))/-8;
    x2(i) = (-38-2*x1(end)+x3(end))/-6;
    x3(i) = (-34+3*x1(end)+x2(end))/7;
    
    E1 = abs((x1(i)-x1(i-1))/x1(i))*100;
    E2 = abs((x2(i)-x2(i-1))/x2(i))*100;
    E3 = abs((x3(i)-x3(i-1))/x3(i))*100;
    
    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002, break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 gauss seidal 법으로 찾은 x1, x2 x3의 값이다')
x1(end)
x2(end)
x3(end)
disp('시행횟수 '), i

% 탐색법
x1 = [0];
x2 = [0];
x3 = [0];

i = 2;
% gauss seidal

while 1
    x1(i) = (-20-x2(i-1)+2*x3(i-1))/-8;
    x2(i) = (-38-2*x1(i-1)+x3(i-1))/-6;
    x3(i) = (-34+3*x1(i-1)+x2(i-1))/7;
    
    E1 = abs((x1(i)-x1(i-1))/x1(i))*100;
    E2 = abs((x2(i)-x2(i-1))/x2(i))*100;
    E3 = abs((x3(i)-x3(i-1))/x3(i))*100;
    
    if E1 < 0.002 && E2 < 0.002 && E3 < 0.002, break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('탐색법 jacobi 법으로 찾은 x1, x2 x3의 값이다')
x1(end)
x2(end)
x3(end)
disp('시행횟수 '), i

% 12-9 비선형 newton raspson 법
i = 1;
xy=[1.2;1.2];


while 1
    J = [2*xy(1,i) 2*xy(2,i);2*xy(1,i) -1];
    f = [xy(1,i)^2+xy(2,i)^2-5; xy(1,i)^2-xy(2,i)+1];
    
    xy(:,i+1) = xy(:,i)- inv(J)*f;
    
    E1 = abs((xy(1,i+1)-xy(1,i))/xy(1,i+1))*100;
    E2 = abs((xy(2,i+1)-xy(2,i))/xy(2,i+1))*100;
    
    if E1 < 0.002 && E2 < 0.002, break
    elseif i>5000, break, end
    
    i = i+1;
end

disp('Newton raspson법을 통하여 얻은값')
x = xy(1,end)
y = xy(2,end)
disp('시행횟수'),i