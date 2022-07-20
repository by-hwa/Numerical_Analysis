clc
clear all
close all

%% 6-20
k1 = 40000;
k2 = 40;
m = 95;
g = 9.81;
h = 0.43;

dd = @(d) (2*k2*d^(5/2))/5 + 0.5*k1*d^2 - m*g*d - m*g*h;

% ezplot(dd,[0,1])
% fzero(dd,0)
% ezplot(fd,[0,1])
% fzero(fd,0)

% 연속 대입법
d(1) =  0;
i = 2;
E = 1;

while E>0.002
    if i>51,break; end
    d(i) = ((2*k2*d(i-1)^(5/2))/(m*g*5)) + ((k1*d(i-1)^2)/(2*m*g)) - h;
    E = ((d(i)-d(i-1))/d(i))*100;
    i = i+1;
end


disp('연속대입법으로 찾은 근은 발산한다')
disp('시도횟수')
disp(length(d))

% 뉴턴법
d = [0];
E = 1;
i = 1;
fd = @(d) k2*d^(3/2) + d*k1 - m*g;

while E>0.002
    if i>51, break, end
    d(i+1) = d(i)- (dd(d(i))/fd(d(i)));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end

disp('뉴턴법으로 찾은 값은')
disp(d(end))
disp('시도횟수')
disp(length(d))

% 할선법
d = [0,0.1];
E = 1;
i = 2;

while E>0.002
    if i>51, break, end
    d(i+1) = d(i) - dd(d(i)) * (d(i-1)-d(i))/(dd(d(i-1))-dd(d(i)));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end

disp('할선법으로 찾은 값은')
disp(d(end))
disp('시도횟수')
disp(length(d))

% 수정할선법

d = [0];
dt = 0.1;
E = 1;
i = 1;

while E>0.002
    if i>51, break, end
    d(i+1) = d(i) - (dd(d(i))*dt)/(dd(d(i)+dt)-dd(d(i)));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end
disp('수정할선법으로 찾은 값은')
disp(d(end))
disp('시도횟수')
disp(length(d))

% brent법
d = [0, 0.1,0.2];
y = [dd(d(1)),dd(d(2)),dd(d(3))];
E = 1;
i = 3;

while E>0.002
    if i>51, break, end
    d(i+1) = ((y(i-1)*y(i))/((y(i-2)-y(i-1))*(y(i-2)-y(i))))*d(i-2) + ((y(i-2)*y(i))/((y(i-1)-y(i-2))*(y(i-1)-y(i))))*d(i-1) + ((y(i-2)*y(i-1))/((y(i)-y(i-2))*(y(i)-y(i-1))))*d(i);
    y(i+1) = dd(d(i+1));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end
disp('brent법으로 찾은 값은')
disp(d(end))
disp('시도횟수')
disp(length(d))


disp('brent<할선법<수정할선법 순으로 시도횟수가 적었다.')
disp('연속대입법은 발산했고 뉴턴법은 값이 이상했다.')


%% 6-21
v0 = 30;
y = 1;
y0 = 1.8;
x = 90;
g = 9.81;
E = 1;

% 연속대입법
th0(1) = 0;
i=2;
while E>0.002
    if i > 51, break,end
    th0(i) = atand((g/(2*v0^2*cosd(th0(i-1))^2))*x - y0/x + y/x);
    E = (th0(i)-th0(i-1))/th0(i)*100;
    i = i+1;
end
disp('6-21번문제')
disp('연속대입법')
disp('초기 세타 값')
disp(th0(end))
disp(length(th0))
disp('번만에 찾았음')

% 뉴턴법
th0 = [0];


 f = @(th0) tand(th0)*x - (g/(2*v0^2*cosd(th0)^2))*x^2 + y0 -y;
df = @(th0) secd(th0)^2*x - (g*sind(th0)/(v0^2*cosd(th0)^3))*x^2;
E =1;
i=1;

while E>0.002
    th0(i+1) = th0(i) - f(th0(i))/df(th0(i));
    E = (th0(i+1)-th0(i))/th0(i+1)*100;
    i = i+1;    
end

disp('')
disp('뉴턴법')
disp('초기 세타 값')
disp(th0(end))
disp(length(th0))
disp('번만에 찾았음')

% 할선법
th0 = [0,0.1];
f = @(th0) tand(th0)*x - (g/(2*v0^2*cosd(th0)^2))*x^2 + y0 -y;
E = 1;
i=2;


while E > 0.002
    if i > 51, break,end
    th0(i+1) = th0(i) - f(th0(i))*((th0(i-1)-th0(i))/(f(th0(i-1)-f(th0(i)))));
    E = (th0(i+1)-th0(i))/th0(i+1)*100;
    i = i+1;
end
disp('')
disp('할선법')
disp('초기 세타 값')
disp(th0(end))
disp(length(th0))
disp('번만에 찾았음')

% 수정할선법

th0 = [0];
f = @(th0) tand(th0)*x - (g/(2*v0^2*cosd(th0)^2))*x^2 + y0 -y;
E = 1;
dth0 = 0.1;
i=1;


while E > 0.002
    if i > 51, break,end
    th0(i+1) = th0(i) - (f(th0(i))*dth0)/(f(th0(i)+dth0)-f(th0(i)));
    E = (th0(i+1)-th0(i))/th0(i+1)*100;
    i = i+1;
end

disp('')
disp('수정할선법')
disp('초기 세타 값')
disp(th0(end))
disp(length(th0))
disp('번만에 찾았음')


 %ezplot(f,[0,90])
% brent법
f = @(th0) tand(th0)*x - (g/(2*v0^2*cosd(th0)^2))*x^2 + y0 ;
th0 = [0,20,40];
y = [f(th0(1)),f(th0(2)),f(th0(3))];
E = 1;
i=3;
gy = @(yr) (((yr-y(i-1))*(yr-y(i)))/((y(i-2)-y(i-1))*(y(i-2)-y(i))))*th0(i-2) + (((yr-y(i-2))*(yr-y(i)))/((y(i-1)-y(i-2))*(y(i-1)-y(i))))*th0(i-1) + (((yr-y(i-2))*(yr-y(i-1)))/((y(i)-y(i-2))*(y(i)-y(i-1))))*th0(i-2);


while E>0.002
    if i > 51, break,end
    th0(i+1) = ((y(i-1)*y(i))/((y(i-2)-y(i-1))*(y(i-2)-y(i))))*th0(i-2) + ((y(i-2)*y(i))/((y(i-1)-y(i-2))*(y(i-1)-y(i))))*th0(i-1) + ((y(i-2)*y(i-1))/((y(i)-y(i-2))*(y(i)-y(i-1))))*th0(i);
    y(i+1) = f(th0(i+1));
    E = abs(((th0(i+1)-th0(i))/th0(i+1))*100);
    i = i+1;
end
disp('brent 법')
disp('초기 세타 값')
disp(th0(end))
disp(length(th0))
disp('번만에 찾았음')

disp('수정할선법<brent법<할선법<연속대입법<뉴턴법 순으로 접근 횟수가 적었다.')
% 6-38
L = 45;
Ec = 1.5*10^11;
Ac = 6.362 * 10^-4;
W = 9000;

% 연속대입법
E = 1;
i = 2;
d(1) = 0;

while E>0.002
    if i > 3000, break,end
    d(i) = ((4500*L)/(Ec*Ac)) + sind(atand(d(i-1)/L))*L;
    E = abs((d(i)-d(i-1))/d(i)) *100;
    i = i+1;
end

disp('6-38문제')
disp('연속대입법')
disp('처짐 d의길이는')
disp(d(end))
disp('시도횟수는')
disp(length(d))
disp('케이블이 늘어난 길이는')
theta = atand(d(end)/L);
delta = (((W/2)/sind(theta))/(Ec*Ac))*L


% 뉴턴법
d = [0];
fd = @(d) ((4500)/(Ec*Ac)) + sind(atand(d/L)) - d/L ;
dfd = @(d) -(2*L*d*cos(L/(d^2+L^2)))/(d^4+2*L^2*d^2+L^4)  - 1/L ;
E = 1;
i = 1;
while E > 0.002
    if i > 100, break,end
    d(i+1) = d(i) - fd(d(i)/dfd(d(i)));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end
disp('뉴턴법')
disp('처짐 d의길이는')
disp(d(end))
disp('시도횟수는')
disp(length(d))
disp('케이블이 늘어난 길이는')
theta = atand(d(end)/L);
delta = (((W/2)/sind(theta))/(Ec*Ac))*L

%할선법
d = [0,0.1];
fd = @(d) ((4500)/(Ec*Ac)) + sind(atand(d/L)) - d/L ;
E = 1;
i=2;

% ezplot(fd,[0,180])
% fzero(fd,1)

while E > 0.002
    if i > 100, break,end
    d(i+1) = d(i) - fd(d(i))*((d(i-1)-d(i))/(fd(d(i-1)-fd(d(i)))));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end

disp('할선법')
disp('처짐 d의길이는')
disp(d(end))
disp('시도횟수는')
disp(length(d))
disp('케이블이 늘어난 길이는')
theta = atand(d(end)/L);
delta = (((W/2)/sind(theta))/(Ec*Ac))*L

% 수정할선법
d = [0];
fd = @(d) ((4500)/(Ec*Ac)) + sind(atand(d/L)) - d/L ;
dt = 0.1;
E = 1;
i = 1;

while E > 0.002
    if i > 51, break,end
    d(i+1) = d(i) - (fd(d(i))*dt)/(fd(d(i)+dt)-fd(d(i)));
    E = abs((d(i+1)-d(i))/d(i+1))*100;
    i = i+1;
end

disp('수정할선법')
disp('처짐 d의길이는')
disp(d(end))
disp('시도횟수는')
disp(length(d))
disp('케이블이 늘어난 길이는')
theta = atand(d(end)/L);
delta = (((W/2)/sind(theta))/(Ec*Ac))*L

% brent 법
fd = @(d) ((4500)/(Ec*Ac)) + sind(atand(d/L)) - d/L ;
d = [0,1,2];
E = 1;
i = 3;
y = [fd(d(1)),fd(d(2)),fd(d(3))];

while E>0.002
    if i > 51, break,end
    d(i+1) = ((y(i-1)*y(i))/((y(i-2)-y(i-1))*(y(i-2)-y(i))))*d(i-2) + ((y(i-2)*y(i))/((y(i-1)-y(i-2))*(y(i-1)-y(i))))*d(i-1) + ((y(i-2)*y(i-1))/((y(i)-y(i-2))*(y(i)-y(i-1))))*d(i);
    y(i+1) = fd(d(i+1));
    E = abs(((d(i+1)-d(i))/d(i+1))*100);
    i = i+1;
end

disp('brent법')
disp('처짐 d의길이는')
disp(d(end))
disp('시도횟수는')
disp(length(d))
disp('케이블이 늘어난 길이는')
theta = atand(d(end)/L);
delta = (((W/2)/sind(theta))/(Ec*Ac))*L
disp('brent법<수정할선법<할선법<연속대입법 순으로 시도횟수가 적었다.')
disp('newton법은 값이 이상하게 나왔다..')

% 연속대입법 뉴턴법 할선법 수정할선법 brent법
% 같은 x0에서 근을 누군 찾고 누군 못찾았는지 누군 몇번만에 찾았는데 다른건 어떤지 비교
