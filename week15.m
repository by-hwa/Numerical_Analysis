clc
clear all
close all
% Ex 23-2

%% ODE

% ydot = @(t,y) 10*exp(-((t-2)^2)/(2*0.075^2)) - 0.6*y;
% ydot(1,1)
% function
%ydot(1,1)

%% Euler
dt = 0.001;
t = [0:dt:4];
L = length(t);

y0 = 0.5;
y(1) = y0;


for i = 1:L-1
    y(i+1) = y(i) + ydot(t(i),y(i))*dt;
end

figure;
plot(t,y, 'b', 'linewidth', 6)
grid on
hold on

dt = 0.1;
t = [0:dt:4];
L = length(t);

y=0;
y0 = 0.5;
y(1) = y0;


for i = 1:L-1
    y(i+1) = y(i) + ydot(t(i),y(i))*dt;
end


plot(t,y, 'r')
grid on

%% Heun

dt = 0.1;
t = [0:dt:4];
L = length(t);

y=0;
y0 = 0.5;
y(1) = y0;
Es = 0.1;

for i = 1:L-1
   yo(i+1) = y(i) + ydot(t(i),y(i))*dt;
   j = 1;
   E = 10;
   yc = yo(i+1);
   while E(j)>Es
       yc(j+1) = y(i) + ((ydot(t(i),y(i)) + ydot(t(i+1),yc(j)))/2)*dt;
       E(j+1) = abs(yc(j+1) - yc(j));
       j = j+1;
   end
   y(i+1) = yc(j);
end

plot(t,y,'k')
grid on

%% RK4
dt = 0.1;
t = [0:dt:4];
L = length(t);

y=0;
y0 = 0.5;
y(1) = y0;

for i = 1:L-1
    k1 = ydot(t(i),y(i));
    k2 = ydot(t(i) + (1/2)*dt, y(i) + (1/2)*k1*dt);
    k3 = ydot(t(i) + (1/2)*dt ,y(i) + (1/2)*k2*dt);
    k4 = ydot(t(i) + dt, y(i)+k3*dt);

    y(i+1) = y(i) + (1/6)*(k1+ 2*k2 + 2*k3 + k4)*dt;
end

plot(t,y,'mo')
grid on

%% ODE 23
dt = 0.1;

y=0;
y0 = 0.5;
E = 0;
y(1) = y0;

DT(1) = dt;
Es = 0.001;
Tor = 0.00000000000001;

T(1) = 0;
i = 1;
while T<=4
    j = 1;
    Next = 0;
    E = 0;
    while Next == 0
        k1 = ydot(T(i),y(i));
        k2 = ydot(T(i) + dt/2, y(i) + k1*dt/2);
        k3 = ydot(T(i) + 3*dt/4, y(i) + 3*k2*dt/4);
        k4 = ydot(T(i) + dt, y(i) + dt);
        y23 = y(i) + (1/9)*(2*k1 + 3*k2 + 4*k3)*dt;
        E(j+1) = abs((1/72)*(-5*k1 + 6*k2 + 8*k3 - 9*k4)*dt);
        if j == 1
            Next = 0;
        else
            if E(j+1) > (Es+Tor) && E(j) > (Es+Tor)
                dt = dt*abs(Es/E(j+1))^0.2;
                Next = 0;
            elseif E(j+1) < (Es-Tor) && E(j) < (Es-Tor)
                dt = dt*abs(Es/E(j+1))^0.25;
                Next = 0;
            else
                Next = 1;
            end
        end
        j = j+1;
    end
    y(i+1) = y23;
    T(i+1) = T(i) + dt;
    DT(i+1) = dt;
    i = i+1;
end

plot(T,y,'go','linewidth',3)

%% ODE 45
% dt = 0.1;
% 
% y=0;
% y0 = 0.5;
% E = 0;
% y(1) = y0;
% 
% DT(1) = dt;
% Es = 0.001;
% Tor = 0.00000000000001;
% 
% T(1) = 0;
% i = 1;
% while T<=4
%     j = 1;
%     Next = 0;
%     E = 0;
%     while Next == 0
%         k1 = ydot(T(i),y(i));
%         k2 = ydot(T(i) + dt/5, y(i) + k1*dt/5);
%         k3 = ydot(T(i) + 3*dt/10, y(i) + 3*k2*dt/40 + 9*k2*dt/40);
%         k4 = ydot(T(i) + dt, y(i) + dt);
%         y45 = y(i) + (1/9)*(2*k1 + 3*k2 + 4*k3)*dt;
%         E(j+1) = abs(y45(i-1)y45(i));
%         if j == 1
%             Next = 0;
%         else
%             if E(j+1) > (Es+Tor) && E(j) > (Es+Tor)
%                 dt = dt*abs(Es/E(j+1))^0.2;
%                 Next = 0;
%             elseif E(j+1) < (Es-Tor) && E(j) < (Es-Tor)
%                 dt = dt*abs(Es/E(j+1))^0.25;
%                 Next = 0;
%             else
%                 Next = 1;
%             end
%         end
%         j = j+1;
%     end
%     y(i+1) = y23;
%     T(i+1) = T(i) + dt;
%     DT(i+1) = dt;
%     i = i+1;
% end

%% MatLAB ODE23
[t,y] = ode23(@ydot, [0,4], y0);

plot(t,y,'ro', 'linewidth',1.5)

figure;
plot(T,DT)
grid on