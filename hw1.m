clc
clear all
close all

disp("문제 8-5")
a = [3+2i 4;-i 1]
b = [2+i;3]

disp("ax = b")
disp("x값을 구하기 위해서 a의 역행렬을 곱해준다")
disp("a-1ax = a-1b")
disp("x = a-1b")

x = inv(a) * b
disp("x값이 올바른지 확인")

a * x
disp("의 값이 b와 같다.")

disp("문제 8-10")
disp("cd = e")
disp("c = Fh1;Fv1;Fh2;Fv2;Fh3;Fv3;")
c = [-cosd(30) 0 cosd(60) 0 0 0;
    -sind(30) 0 -sind(60) 0 0 0;
    cosd(30) 1 0 0 0 1;
    sind(30) 0 0 1 0 0;
    0 -1 -cosd(60) 0 0 0;
    0 0 sind(60) 0 1 0]

e = [0;1000;0;0;0;0;]

disp("c-1cd = c-1e")
disp("d = c-1e")
disp("d = F1;F2;F3;V2;V3;H2")

d = inv(c) * e


disp("문제 8-11")
disp("k = 계수행렬 x = x 행렬 a = 가속도 행렬")

syms k1 k2 k3 k4 m1 m2 m3 x1 x2 x3
k = [-(k1+k2)/m1 k2/m1 0;k2/m2 -(k2+k3)/m2 k3/m2;0 k3/m3 -(k3+k3)/m3]

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

a = k * x


disp("문제 8-13")
g = 9.81;
m1 = 2;
m2 = 3;
m3 = 2.5;
k = 10;

disp("K행렬 * 변위행렬 = mg")
ka = [2*k -k 0;k 0 -k;0 k -k]
m = [m1*g;m2*g;m3*9]

disp("ka x = m")
disp("ka-1ka x = ka-1m")

x = inv(ka) * m


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

disp("c i = v")
disp("c-1ci=c-1v")
disp("i=c-1v")

i = inv(c) * v

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

i = inv(c) * v