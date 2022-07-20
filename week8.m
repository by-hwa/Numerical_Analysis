clc
clear all
close all

%% ex1
Fx = @(x) 1*x^2 - 2*x+1;

% ezplot(Fx)

fzero(Fx,2)

%% fixed- point iteration
x0 = 0;
Es = 0.0000001;

Ea = 1;

x = x0;
i = 1;

while Ea(i) > Es
    x(i+1) = 0.5 * x(i)^2 + 0.5;
    Ea(i+1) = abs(x(i+1)-x(i)); 
    i = i+1;
end

i
x_fixed = x(i)
Ea_fixed = Ea(i)

figure
plot(x)

%% Newton-Raphson
x0 = 0.5;
Es = 0.01;
Ea = 1;

x = x0;
Fdotx = @(x) 2*x-2;

f = Fx(x);
fdot = Fdotx(x);
i = 1

while Ea(i) > Es
    x(i+1) = x(i) - f(i)/fdot(i);
    f(i+1) = Fx(x(i+1));
    fdot(i+1) = Fdotx(x(i+1));
    Ea(i+1) = abs(x(i+1)-x(i)); 
    i = i+1;
end

i
x_fixed = x(i)
Ea_fixed = Ea(i)

%% secent

x0 = 0.5;
Es = 0.01;
Ea = 1;

x = x0;
f = Fx(x);

i = 1;
while Ea(i) > Es
    x(i+1) = x(i) - (f(i));
    f(i+1) = Fx(x(i+1));
    Ea(i+1) = abs(x(i+1)-x(i)); 
    i = i+1;
end

%% inverse 2nd interpolation(역2차 보간법 brent법)

