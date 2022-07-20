clc
clear all
close all
%% 22.2
fdot = @(x,y) (1 + 2 * x) * sqrt(y);
dt = 0.25;
X = [0:dt:1];
y0 = 1;

% 해석적 기법
% Euler 법
fdot = @(x,y) (1 + 2 * x) * sqrt(y);
dt = 0.25;
X = [0:dt:1];
L = length(X);
y0 = 1;
y(1) = y0;
for i =1:L-1
    y(i+1) = y(i) + fdot(X(i),y(i))*dt;
end

figure;
plot(X,y)
hold on

% 반복 없는 Heun 법
dt = 0.25;
X = [0:dt:1];
L = length(X);
y=0;

y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(X(i),y(i));
    k2 = fdot(X(i)+dt,y(i)+k1*dt);
    y(i+1) = y(i) + ((k1/2)+(k2/2))*dt;
end

plot(X,y)

% Ralston 법
dt = 0.25;
X = [0:dt:1];
L = length(X);
y=0;

y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(X(i),y(i));
    k2 = fdot(X(i)+dt*2/3,y(i)+k1*dt*2/3);
    y(i+1) = y(i) + ((k1/4)+(k2*3/4))*dt;
end

plot(X,y)

% 4차 RK법
dt = 0.25;
X = [0:dt:1];
L = length(X);
y=0;

y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(X(i),y(i));
    k2 = fdot(X(i)+dt/2,y(i)+k1*dt/2);
    k3 = fdot(X(i)+dt/2,y(i)+k2*dt/2);
    k4 = fdot(X(i)+dt,y(i)+k3*dt);
    y(i+1) = y(i) + (1/6)*(k1+ 2*k2 + 2*k3 + k4)*dt;
end

plot(X,y)
legend({'Euler','Heun','Ralston','4차 RK'})
%% 22.3

fdot = @(t,y) -y + t^2;

% 수정반복 없는 Heun법
dt = 0.5;
t = [0:dt:3];
L = length(t);
y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(t(i),y(i));
    k2 = fdot(t(i)+dt,y(i)+k1*dt);
    y(i+1) = y(i) + ((k1/2)+(k2/2))*dt;
end

figure;
plot(t,y);
hold on

% 수정반복 있는 Heun 법 Es = 0.1%
dt = 0.5;
t = [0:dt:3];
L = length(t);
y0 = 1;
y(1) = y0;
Es = 0.1;

for i = 1:L-1
    j = 1;
    yo = y(i) + fdot(t(i),y(i))*dt;
    E = 1;
    yc(j) = yo;
    while E>Es
        yc(j+1) = y(i) + ((fdot(t(i),y(i)) + fdot(t(i+1),yc(j)))/2)*dt;
        E = abs((yc(j)-yc(j+1))/yc(j))*100;
        j = j+1;
    end
    y(i+1) = yc(j);
end

plot(t,y);

% 중점법
dt = 0.5;
t = [0:dt:3];
L = length(t);
y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(t(i),y(i));
    k2 = fdot(t(i)+dt/2,y(i)+k1*dt/2);
    y(i+1) = y(i) + k2*dt;
end


plot(t,y);

% Ralston 법

dt = 0.5;
t = [0:dt:3];
L = length(t);
y0 = 1;
y(1) = y0;

for i = 1:L-1
    k1 = fdot(t(i),y(i));
    k2 = fdot(t(i)+dt*2/3,y(i)+k1*dt*2/3);
    y(i+1) = y(i) + ((k1/4)+(k2*3/4))*dt;
end

plot(t,y);
legend({'반복없는 Heun','반복있는 Heun', '중점법', 'Ralston법'})
%% 22.9
dt = 0.1;
t = [0:dt:4];
L = length(t);


y0 = 1;
ydot0 = 0;

yddot = @(t, y, yd) -9*y;
z = @(t, y, yd) yd;

y(1) = y0;
yd(1) = ydot0;

% Euler
for i = 1:L-1
    yd(i+1) = yd(i) + yddot(t(i), y(i), yd(i))*dt;
    y(i+1) = y(i) + yd(i+1)*dt;
end

figure;
hold on
f = @(t) cos(3*t);
ezplot(f,[0,4]);
plot(t,y);


% 4차 RK 법
y = 0;
yd = 0;

y0 = 1;
ydot0 = 0;

y(1) = y0;
yd(1) = ydot0;

for i = 1:L-1
    k11 = z(t(i),y(i),yd(i));
    k12 = yddot(t(i),y(i),yd(i));
    
    k21 = z(t(i)+dt/2,y(i)+k11*dt/2,yd(i)+k12*dt/2);
    k22 = yddot(t(i)+dt/2,y(i)+k11*dt/2,yd(i)+k12*dt/2);
    
    k31 = z(t(i)+dt/2,y(i)+k21*dt/2,yd(i)+k22*dt/2);
    k32 = yddot(t(i)+dt/2,y(i)+k21*dt/2,yd(i)+k22*dt/2);
    
    k41 = z(t(i)+dt,y(i)+k31*dt,yd(i)+k32*dt);
    k42 = yddot(t(i)+dt,y(i)+k31*dt,yd(i)+k32*dt);
    
    y(i+1) = y(i) + (1/6)*(k11+ 2*k21 + 2*k31 + k41)*dt;
    yd(i+1) = yd(i) + (1/6)*(k12+ 2*k22 + 2*k32 + k42)*dt;
end

plot(t,y)

% [t,y] = ode23(fydot, [0,4], y0);
% 
% plot(t,y,'ro', 'linewidth',1.5)