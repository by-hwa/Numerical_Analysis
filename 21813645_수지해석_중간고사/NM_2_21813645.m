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
        
    if j > 20, break, end
    
    x(i+2) = x(i) + dt;
    i = i+1;
end

dd = [];

j = 1;
for i = 1:20
    d1 = d(i);
    d2 = d(i+1);
    while 1
        dis = (1.6180339-1)*(d2-d1);
        dx1 = d1+dis;
        dx2 = d2-dis;
        
        E = 0.3820 * abs(d2-d1);

        if E < 0.000001
            if f(dx2)<0
                dd(j)= dx2; j=j+1;
            end
            break
        end
        
        if f(dx1)>f(dx2)
            d1 = dx2;
        elseif f(dx1)<f(dx2)
            d2 = dx1;
        end
    end
end

left = [];
right = [];
j = 1;

for i = 1:2:10
    right(j) = dd(i);
    j = j+1;
end

j = 1;
for i = 2:2:10
    left(j) = dd(i);
    j = j+1;
end

disp('시작점이')
now
disp('일때 최적점 위치')
left
right

