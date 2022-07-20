clc
clear all
close all

%%

m = 1;
c = 0.3;
c1 = 56.7887357729738;
c2 = -49.783067365588;
c3 = -0.474875386694339;
c4 = -0.00566840738588095;
w = 25.1327412287183;

x0 = 7;
v0 = 3;

dt = 0.05;

T = 2;
t = [0:dt:T];
L = length(t);
F = 300*sin(w*t);

fx = @(t) c1 + c2 *exp(-c/m*t) + c3*sin(w*t) + c4*cos(w*t);

% figure;
% plot(t,F)

%% 1. 질량의 위치 2초를 0.05 간격으로 그래프.

for i = 1:L
    x(i) = fx(t(i));
end

figure;
plot(t,x)
title('1. 질량위치 그래프')
%% 2. 위치데이터를 바탕으로 질량의 속도를 계산 후, 그래프

v(1) = v0;
for i = 2:L
    v(i) = (fx(i) - fx(i-1))/dt;
end

figure;
plot(t,v)
title('2. 속도 계산 그래프')
%% 3. 운동방정식을 바탕으로 속도 계산, 그래프

ddx = @(t,x,v) (300*sin(w*t)-(c*v))/m;
dx = @(t,x,v) v;

% 고정식
x(1) = x0;
v(1) = v0;

for i = 1:L-1
    
    k11 = dx(t(i),x(i),v(i));
    k12 = ddx(t(i),x(i),v(i));
    
    k21 = dx(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2);
    k22 = ddx(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2);
    
    k31 = dx(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2);
    k32 = ddx(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2);
    
    k41 = dx(t(i)+dt,x(i)+(k31*dt),v(i)+k32*dt);
    k42 = ddx(t(i)+dt,x(i)+(k31*dt),v(i)+k32*dt);
    
    wx(i+1) = x(i) + (1/6)*(k11+2*k21+2*k31+k41);
    wv(i+1) = v(i) + (1/6)*(k12+2*k22+2*k32+k42);
    
end

figure;
plot(t,wx);
title('3.운동방정식 위치 그래프')
figure;
plot(t,wv)
title('3.운동 방정식 속도 그래프')

% 적응식
% 
% tt(i) = 0;
% TT(i) = 0;
% Es = 0.0000001;
% j = 1;
% E = 1;
% h = 0.05;
% while TT(end) < 2
%     
%     while E(end)>Es
%        k11 = dx(t(i),x(i),v(i));
%        k12 = ddx(t(i),x(i),v(i));
%     
%        k21 = dx(t(i)+h/2,x(i)+k11*h/2,v(i)+k12*h/2);
%        k22 = ddx(t(i)+h/2,x(i)+k11*h/2,v(i)+k12*h/2);
%     
%        k31 = dx(t(i)+h*3/4,x(i)+k21*h*3/4,v(i)+k22*h*3/4);
%        k32 = ddx(t(i)+h*3/4,x(i)+k21*h*3/4,v(i)+k22*h*3/4);
%     
%        awx(i+1) = x(i) + (1/9)*(2*k11+3*k21+4*k31);
%        awv(i+1) = v(i) + (1/9)*(2*k12+3*k22+4*k32);
%     
%        E(i+1) = (1/72)*(-5*k12+6*k22+8*k32-9*k42);
%        
%        if E(i+1)>Es && E(i)>Es
%            h = h*(Es/E(i+1))^0.2;
%        elseif E(i+1)<Es && E(i)<Es
%            h = h*(Es/E(i+1))^0.25;
%        else
%            break
%        end
%     end
%     TT(i+1) = TT(i) + h;
% end

    
        
    

%% 4. 그래프 비교하기

% 위치
figure;
plot(t,wx)
title('4.위치 그래프 비교')
hold on
plot(t,x)
legend({'운동방정식','데이터 점'})
% 속도
figure;
plot(t,wv)
title('4.속도 그래프 비교')
hold on
plot(t,v)
legend({'운동방정식','데이터 점'})

%% 5.속도가 양수에서 음수로 바뀌는 것을 감지하면
% 고정식
x(1) = x0;
v(1) = v0;

ddx = @(t,x,v,c) (300*sin(w*t)-(c*v))/m;
dx = @(t,x,v) v;

c(1) = 0.3;
j = 2;

for i = 1:L-1
    
    k11 = dx(t(i),x(i),v(i));
    k12 = ddx(t(i),x(i),v(i),c(end));
    
    k21 = dx(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2);
    k22 = ddx(t(i)+dt/2,x(i)+k11*dt/2,v(i)+k12*dt/2,c(end));
    
    k31 = dx(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2);
    k32 = ddx(t(i)+dt/2,x(i)+k21*dt/2,v(i)+k22*dt/2,c(end));
    
    k41 = dx(t(i)+dt,x(i)+(k31*dt),v(i)+k32*dt);
    k42 = ddx(t(i)+dt,x(i)+(k31*dt),v(i)+k32*dt,c(end));
    
    wx5(i+1) = x(i) + (1/6)*(k11+2*k21+2*k31+k41);
    wv5(i+1) = v(i) + (1/6)*(k12+2*k22+2*k32+k42);
    
    if wv5(i) * wv5(i+1)<0
        if wv5(i) > 0
            c(j) = c(j-1)/2;
            j = j+1;
        end
    end
    
end

figure;
plot(t,wv5)
title('5.제어된 속도 그래프')

figure;
plot(c)
title('5.변경된 뎀퍼 상수 그래프')