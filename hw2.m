clc
clear all
close all

disp("문제 8-5")
A = [3+2i 4;-i 1]
B = [2+i;3]

disp('Ax = B')
disp('C = [A B]')

C = [A B]

[m,n] = size(A);

disp('1) 가우스 소거 + 후방대입')
disp('gauss 소거')
for i = 1:m
    for j=(i+1):m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i))*C(i,:)
    end
end


x = zeros(n,1);

disp('x에 후방대입')
for i = n:-1:1
    x(i,:) = (C(i,n+1) - C(i,1:n)*x)/C(i,i)
end

disp('A=LU, LUx = B, Ux = Y, LY = B')
disp('LY = B 전방대입, Ux=Y 후방대입')
disp('2) LU 분해 후 전방대입 -> 후방대입')

[L,U] = lu(A)

disp('Y 전방대입')

Y = zeros(n,1);

for i = 1:n
    Y(i,:) = (B(i,:) - L(i,1:n)*Y)/L(i,i)
end

disp('U 후방대입')

x = zeros(n,1);

for i = n:-1:1
    x(i,:) = (Y(i,:) - U(i,1:n)*x)/U(i,i)
end

disp("문제 8-10")
disp("Ax = B")
disp("A = Fh1;Fv1;Fh2;Fv2;Fh3;Fv3;")
A = [-cosd(30) 0 cosd(60) 0 0 0;
    -sind(30) 0 -sind(60) 0 0 0;
    cosd(30) 1 0 0 0 1;
    sind(30) 0 0 1 0 0;
    0 -1 -cosd(60) 0 0 0;
    0 0 sind(60) 0 1 0]

B = [0;1000;0;0;0;0;]

C = [A B]


[m,n] = size(A);
disp('1) 가우스 소거 + 후방대입')
disp('gauss 소거를 할 수있도록 행렬을 재배치 해준다.')
C = [C(1,:);C(5,:);C(2,:);C(4,:);C(6,:);C(3,:)]
A = [C(:,1:6)]
B = [C(:,7)]

disp('gauss 소거')

for i = 1:m
    for j = (i+1):m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i)*C(i,:))
    end
end

x = zeros(n,1);
disp('후방대입')

for i = n:-1:1
    x(i,:) = (C(i,n+1)-C(i,1:n)*x)/C(i,i)
end

disp('2) LU 분해 후 전방대입 -> 후방대입')
disp('LY = B 전방대입, Ux=Y 후방대입')
[L,U] = lu(A)

disp('Y 전방대입')

Y = zeros(n,1);

for i = 1:n
    Y(i,:) = (B(i,:) - L(i,1:n)*Y)/L(i,i)
end

disp('U 후방대입')

x = zeros(n,1);

for i = n:-1:1
    x(i,:) = (Y(i,:) - U(i,1:n)*x)/U(i,i)
end

disp("문제 8-11")
disp("k = 계수행렬 x = x 행렬 a = 가속도 행렬")


k1 = 10;
k2 = 30;
k3 = 30;
k4 = 10;
m1 = 1;
m2 = 1;
m3 = 1;

k = [-(k1+k2)/m1 k2/m1 0;k2/m2 -(k2+k3)/m2 k3/m2;0 k3/m3 -(k3+k3)/m3]

x1 = 0.05;
x2 = 0.04;
x3 = 0.03;

x = [x1;x2;x3]

disp('inv 연산이 필요없다')
a = k*x

disp("문제 8-13")
g = 9.81;
m1 = 2;
m2 = 3;
m3 = 2.5;
k = 10;

disp("K행렬 * 변위 = mg")
ka = [2*k -k 0;k 0 -k;0 k -k]
mg = [m1*g;m2*g;m3*9]

disp('하삼각과 상삼각을 만들기 위해 2행과 3행을 교환')
ka = [ka(1,:);ka(3,:);ka(2,:)]
mg = [mg(1,:);mg(3,:);mg(2,:)]

[m,n] = size(ka);

C = [ka mg]

disp('1) 가우스 소거 + 후방대입')
disp('gauss 소거')


for i=1:m
    for j =(i+1):m
        C(j,:) = C(j,:)-(C(j,i)/C(i,i))*C(i,:)
    end
end

x = zeros(n,1);
disp('후방대입')

for i = n:-1:1
    x(i,:) = (C(i,n+1)-C(i,1:n)*x)/C(i,i)
end

disp('2) LU 분해 후 전방대입 -> 후방대입')
disp('ka=LU, LUx = mg, Ux = Y, LY = mg')
disp('LY = mg 전방대입, Ux=Y 후방대입')

[L, U] = lu(ka)

disp('Y전방대입')
Y = zeros(n,1);

for i = 1:n
    Y(i,:) = (mg(i,:) - L(i,1:n)*Y)/L(i,i)
end

x = zeros(n,1);

disp('U후방대입')
for i = n:-1:1
    x(i,:) = (Y(i,:)-U(i,1:n)*x)/U(i,i)
end

disp("문제 8-14")
disp("i12-i23-i25 = 0")
disp("i45-i56+i25 = 0")
disp("i23-i34 = 0")
disp("i34-i45 = 0")
disp("r23i23+r34i34+r45i45-r25i25=0")
disp("r12i12+r25i25+r56i56=150")


disp("i = [i12;i23;i34;i45;i56;i25]")

r12 = 5;
r23 = 20;
r34 = 2;
r45 = 5;
r56 = 25;
r25 = 5;

c = [1 -1 0 0 0 -1;
    0 0 0 1 -1 1;
    0 1 -1 0 0 0;
    0 0 1 -1 0 0;
    0 r23 r34 r45 0 -r25;
    r12 0 0 0 r56 r25;]
v = [0;0;0;0;0;150]

disp('LU를 상삼각 하삼각으로 만들기 위해 이동')
c = [c(6,:);c(5,:);c(3,:);c(4,:);c(1,:);c(2,:)]
v = [v(6,:);v(5,:);v(3,:);v(4,:);v(1,:);v(2,:)]


[m,n] = size(c)

C = [c v]

disp('1) 가우스 소거 + 후방대입')
disp('gauss 소거')

for i = 1 : m
    for j = i+1:m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i))*C(i,:)
    end
end

i = zeros(n,1);
disp('후방대입')

for j = n:-1:1
    i(j,:) = (C(j,n+1) - (C(j,1:n)*i))/C(j,j)
end

disp('2) LU 분해 후 전방대입 -> 후방대입')
disp('c=LU, LUx = v, Ux = Y, LY = v')
disp('LY = v 전방대입, Ux=Y 후방대입')

[L, U] = lu(c)

disp('Y전방대입')
Y = zeros(n,1);

for i = 1:n
    Y(i,:) = (v(i,:) - L(i,1:n)*Y)/L(i,i)
end

x = zeros(n,1);

disp('U후방대입')
for i = n:-1:1
    x(i,:) = (Y(i,:)-U(i,1:n)*x)/U(i,i)
end

x = zeros(n,1);

disp('U후방대입')
for i = n:-1:1
    x(i,:) = (Y(i,:)-U(i,1:n)*x)/U(i,i)
end

disp("문제 8-15")
disp("node2 i12-i23-i25 = 0")
disp("node3 i23-i34-i35 = 0")
disp("node4 i34 - i45 = 0")
disp("node5 i45-i56+i25+i35 = 0")
disp("loop1 r12i12+r56i56+r25i25 = -120")
disp("loop2 r23i23+r35i35-r25i25 = 0")
disp("loop3 r34i34+r45i45-r35i35 = 0")
disp("i = [i12;i23;i34;i45;i56;i25;i35]")

r12 = 35;
r23 = 30;
r34 = 8;
r45 = 15;
r56 = 5;
r25 = 10;
r35 = 7;

c = [1 -1 0 0 0 -1 0;
    0 1 -1 0 0 0 -1;
    0 0 1 -1 0 0 0;
    0 0 0 1 -1 1 1;
    r12 0 0 0 r56 r25 0;
    0 r23 0 0 0 -r25 r35;
    0 0 r34 r45 0 0 -r35]

v = [0;0;0;0;-120;0;0]

c = [c(5,:);c(6,:);c(7,:);c(3,:);c(4,:);c(1,:);c(2,:)]
v = [v(5,:);v(6,:);v(7,:);v(3,:);v(4,:);v(1,:);v(2,:)]

[m,n] = size(c);



C = [c v]
disp('1) 가우스 소거 + 후방대입')
disp('gauss 소거')

for i = 1 : m
    for j = i+1:m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i))*C(i,:)
    end
end

i = zeros(n,1);
disp('후방대입')

for j = n:-1:1
    i(j,:) = (C(j,n+1) - (C(j,1:n)*i))/C(j,j)
end

[L, U] = lu(c)

disp('Y전방대입')
Y = zeros(n,1);

for i = 1:n
    Y(i,:) = (v(i,:) - L(i,1:n)*Y)/L(i,i)
end

x = zeros(n,1);

disp('U후방대입')
for i = n:-1:1
    x(i,:) = (Y(i,:)-U(i,1:n)*x)/U(i,i)
end

x = zeros(n,1);

disp('U후방대입')
for i = n:-1:1
    x(i,:) = (Y(i,:)-U(i,1:n)*x)/U(i,i)
end