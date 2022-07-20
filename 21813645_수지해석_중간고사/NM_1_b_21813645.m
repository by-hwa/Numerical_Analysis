clc
clear all
close all

f  = @(x) (1/1.6499999952204)*(sin(0.3*pi*(x+5.91589626160177))+0.65);

d = [];

dt = 0.1;

now = 0; % 시작점

x = [now,dt];
i = 1;
j = 1;


while 1
    f1 = f(x(i));
    f2 = f(x(i+1));
    
    if f1 * f2 <0
        k = 2;
        xx = [(x(i)+x(i+1))/2];
        while 1 
            xr = x(i+1) - (f(x(i+1))*(x(i)-x(i+1)))/(f(x(i))-f(x(i+1)));
            xx(k) = xr;
            E = abs((xx(k-1)-xx(k))/xx(k-1))*100;
            k = k+1;
            if E<0.000001
                j = j + 1;
                d(j) = xx(k-1);
                break
            end
        end
    end
        
    if j > 10, break, end
    
    x(i+2) = x(i) + dt;
    i = i+1;
end

left = [];
right = [];
j = 1;

for i = 3:2:11
    right(j) = d(i);
    j = j+1;
end

j = 1;
for i = 2:2:10
    left(j) = d(i);
    j = j+1;
end

disp('시작점이')
now
disp('일때 근의 위치')
left
right
