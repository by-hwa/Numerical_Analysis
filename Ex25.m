clc
clear all
close all
%% Ex 22.5

cd = 0.25;
g = 9.81;
m = 68.1;

xd = @(t,x,v) v;
vd = @(t,x,v) g-cd/m*v^2;

dt = 2;
t = [0:dt:10];
L = length(t);

x0 = 0;
v0 = 0;

x(1) = x0;
v(1) = v0;

for i = 1:L-1
    k11 = xd(t(i),x(i),v(i));
    k12 = vd(t(i),x(i),v(i));

    k21 = xd(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2);
    k22 = vd(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2);
    
    k31 = xd(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2);
    k32 = vd(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2);
    
    k41 = xd(t(i)+dt,x(i)+k31*dt,v(i)+k31*dt);
    k42 = vd(t(i)+dt,x(i)+k31*dt,v(i)+k31*dt);
    
    
    break
end