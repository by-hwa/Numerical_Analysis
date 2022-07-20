clc
clear all
close all
%% 5-11
N = [0:0.1:12];

for i = 1:length(N)
    if N(i)<3
        M(i) = - (100/18)*N(i)^3 + 265*N(i);
    elseif N(i)>=3 && N(i)<6
        M(i) = - 50*N(i)^2 + 415*N(i) - 150;
    elseif N(i)>=6 && N(i)<10
        M(i) = - 185*N(i) + 1650;
    else
        M(i) = 100*N(i) - 1200;
    end
end
figure(Name = 'Moment 선도')
plot(N,M)

% 증분법
xL = 0;
xU = 12;
N = 50;

dx = (xU-xL)/N;
x = [xL:dx:xU];
j = 1

for i = 1:length(x)-1
    if x(i)<3
        M(i) = - (100/18)*x(i)^3 + 265*x(i);
        M(i+1) = - (100/18)*x(i+1)^3 + 265*x(i+1);
    elseif x(i)>=3 && x(i)<6
        M(i) = - 50*x(i)^2 + 415*x(i) - 150;
        M(i+1) = - 50*x(i+1)^2 + 415*x(i+1) - 150;
    elseif x(i)>=6 && x(i)<10
        M(i) = - 185*x(i) + 1650;
        M(i+1) = - 185*x(i+1) + 1650;
    else
        M(i) = 100*x(i) - 1200;
        M(i+1) = 100*x(i+1) - 1200;
    end
    S = M(i) * M(i+1);
    if S == 0
            if M(i) == 0
            xr(j) = x(i);
            j = j+1;
        else
            xr(j) = x(i+1);
            j=j+1;
        end
    elseif S<0
        xr(j) = abs(x(i) + x(i+1))/2;
        j = j+1;
    end
end

disp('근을 갖는 값은')
xr

% 함수 선언
mf1 = @(x) - (100/18)*x^3 + 265*x;
mf2 = @(x) - 50*x^2 + 415*x - 150;
mf3 = @(x)  - 185*x + 1650;
mf4 = @(x) 100*x - 1200;

% 이분법
xL = 0;
xU = 12;
Es = 0.001;
N = log2((xU-xL)/Es);
N = ceil(N);
xr = [0];

for i = 1:N
    xr(i) = (xU+xL)/2;
    if xL<3, ML = mf1(xL);
    elseif xL>=3 && xL<6, ML = mf2(xL);
    elseif xL>=6 && xL<10, ML = mf3(xL);
    else ML = mf4(xL);end
    if xr(i)<3, MR = mf1(xr(i));
    elseif xr(i)>=3 && xr(i)<6, MR = mf2(xr(i));
    elseif xr(i)>=6 && xr(i)<10, MR = mf3(xr(i));
    else MR = mf4(xr);end
    if xU<3, MU = mf1(xU);
    elseif xU>=3 && xU<6, ML = mf2(xU);
    elseif xU>=6 && xU<10, ML = mf3(xU);
    else MU = mf4(xU);end
    
    SL = ML * MR;
    SU = MU * MR;
    

    if SL<0
        xU = xr(i);
    else
        xL = xr(i);
    end
end

disp('이분법으로 찾은 근은')
xr(1,end)
disp('이다')

%선형 보간법
xL = 0;
xU = 12;
Es = 0.001;
N = log2((xU-xL)/Es);
N = ceil(N);
xr = [0];
for i = 1:N
    
    if xU<3, MU = mf1(xU);fu=MU;
    elseif xU>=3 && xU<6, MU = mf2(xU);fu=MU;
    elseif xU>=6 && xU<10, MU = mf3(xU);fu=MU;
    else MU = mf4(xU);fu=MU;end
    if xL<3, ML = mf1(xL);fl=ML;
    elseif xL>=3 && xL<6, ML = mf2(xL);fl=ML;
    elseif xL>=6 && xL<10, ML = mf3(xL);fl=ML;
    else ML = mf4(xL);fl=ML;end
    
    if fl == 0 || fu == 0
        xr(i) = (xU+xL)/2;
    else
        xr(i) = xU - (fu*(xU-xL))/(fl-fu);
    end

    if xr(i)<3, MR = mf1(xr(i));
    elseif xr(i)>=3 && xr(i)<6, MR = mf2(xr(i));
    elseif xr(i)>=6 && xr(i)<10, MR = mf3(xr(i));
    else MR = mf4(xr(i));end


    SL = ML * MR;
    SU = MU * MR;
    
    
    if SL<0
        xU = xr(i);
    else
        xL = xr(i);
    end
end

disp('선형보간법으로 찾은 근은')
xr(1,end)
disp('이다')

%% 5-15
disp('5-15')
% 함수 선언
L = 600;
E = 5000;
I = 30000;
w0 = 2.5*10^3;

y = @(x) (w0/(120*E*I*L))*(-x^5+2*L^2*x^3-L^4*x);

figure;
ezplot(y,[0,600])

% 증분법
xL = 0;
xU = L;

N = 500;
dt = (xU-xL)/N;

x = [xL:dt:xU];
j=1;

for i = 1:length(x)
    Y(i) = y(x(i));
    Y(i+1) = y(x(i));
 
    S = Y(i) * Y(i+1);
    
    if S == 0
        if Y(i) == 0
            xr2(j) = x(i);
            j = j+1;
        else
            xr2(j) = x(i+1);
            j=j+1;
        end
    elseif S<0
        xr2(j) = abs(x(i)+x(i+1))/2;
        j = j+1;
    end
end
disp('증분법으로 찾은 근은')
xr2

%이분법

xL = 0;
xU = L;
Es = 0.001;
N = log2((xU-xL)/Es);
N = ceil(N);
xr2 = [0];

for i = 1:N
    xr2(i) = (xU+xL)/2;
    
    EL = y(xL);
    ER = y(xr2(i));
    EU = y(xU);
    
    SL2 = EL * ER;
    SU2 = EU * ER;
    
    
    if SL2<0
        xU = xr2(i);
    else
        xL = xr2(i);
    end
end

disp('이분법으로 찾은 근은')
xr2(1,end)
disp('이다')
disp('0~600사이의 근이 없어서 제일 끝의 600근처의 근을 찾았다.')

%선형 보간법
xL = 0;
xU = L;
Es = 0.001;
N = log2((xU-xL)/Es);
N = ceil(N);
xr2 = [0];
for i = 1:N
    fel = y(xL);
    feu = y(xU);
    
    if fel == 0 || feu == 0
        xr2(i) = (xU+xL)/2;
    else
        xr2(i) = xU - (feu*(xU-xL))/(fel-feu);
    end
    
    EL = y(xL);
    ER = y(xr2(i));
    EU = y(xU);

    SL2 = EL * ER;
    SU2 = EU * ER;
    
    
    if SL2<0
        xU = xr2(i);
    else
        xL = xr2(i);
    end
end

disp('선형보간법으로 찾은 근은')
xr2(1,end)
disp('이다')
%% 5-22

u = 1800;
m0 = 160000;
q = 2600;
g=9.81;

v = @(t) u*log(m0/(m0-q*t)) - g*t
figure;
ezplot(v,[10,50])

% 증분법

ti = 10;
tf = 50;
N = 50;
dt = (tf-ti)/N;
xt = [ti:dt:tf];

j=1;

for i = 1:length(xt)-1
    V(i) = v(xt(i)) - 750;
    V(i+1) = v(xt(i+1)) - 750;
    
    S = V(i) * V(i+1);
    
    if S == 0;
        if V(i) == 0
            xtr(j) = xt(i);
            j = j+1;
        else
            xtr(j) = xt(i+1);
            j = j +1;
        end
    elseif S<0
        xtr(j) = abs(xt(i)+xt(i+1))/2;
        j = j+1;
    end
end

disp('증분법으로 구한 v가 750m/s에 도달하는 시간은')
xtr

% 이분법

ti = 10;
tf = 50;
Es = 0.001;
N = log2((tf-ti)/Es);
N = ceil(N);
xtr = [0];

for i = 1:N
    xtr(i) = (tf+ti)/2;
    
    VL = v(ti)-750;
    VR = v(xtr(i))-750;
    VU = v(tf)-750;
    
    SL3 = VL * VR;
    SU3 = VU * VR;
    
    
    if SL3<0
        tf = xtr(i);
    else
        ti = xtr(i);
    end
end

disp('이분법으로 구한 v가 750m/s에 도달하는 시간은')
xtr(1,end)
        
% 선형보간법

ti = 10;
tf = 50;
Es = 0.001;
N = log2((tf-ti)/Es);
N = ceil(N);
xtr = [0];

for i = 1:N
    
    fi = v(ti)-750;
    ff = v(tf)-750;
    
    if fi == 0 || ff == 0
        xtr(i) = (ti+tf)/2;
    else
        xtr(i) = tf - (ff*(ti-tf))/(fi-ff);
    end

    
    VL = v(ti)-750;
    VR = v(xtr(i))-750;
    VU = v(tf)-750;
    
    SL3 = VL * VR;
    SU3 = VU * VR;
    
    
    if SL3<0
        tf = xtr(i);
    else
        ti = xtr(i);
    end
end

disp('선형보간법으로 구한 v가 750m/s에 도달하는 시간은')
xtr(1,end)