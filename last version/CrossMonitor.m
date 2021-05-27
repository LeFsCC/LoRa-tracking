% x = -2:0.05:2;
% y = -0.1:-0.02:-2.1;

% x = [-4:0.1:-2 ones(1, 20)*-2 -2:0.1:2 ones(1, 10)*2 2:-0.1:1];
% y = [-0.5:-0.1:-4.5 -4.5:0.1:-0.5 -0.5:-0.1:-3.5];

% x = [-4:0.1:-2 ones(1, 20)*-2 -2:-0.1:-4 ones(1, 20)*-4];
% y = [ones(1, 20)*-2 -2:-0.1:-4 ones(1, 20)*-4 -4:0.1:-2];

% x = [0.1:0.1:2.1 ones(1, 20)*2.1 2.1:-0.1:0.1 ones(1, 20)*0.1];
% y = [ones(1, 20)*-0.1 -0.1:-0.1:-2.1 ones(1, 20)*-2.1 -2.1:0.1:-0.1];

% x = -1.0:0.02:1.0;
% y = ones(1, 101)*-1.0;
% 
% x = [-1.0:0.04:0 0:-0.05:-1.0];
% y = [-1.0:-0.02:-1.50 -1.50:-0.02:-2.0];

% x = [-10:0.1:-8 ones(1, 20)*-8 -8:-0.1:-10 ones(1, 20)*-10];
% y = [ones(1, 20)*-10 -10:-0.1:-12 ones(1, 20)*-12 -12:0.1:-10];

% x = [-2.0:0.1:-0.1 ones(1, 20)*-0.1 -0.1:-0.1:-2.0 ones(1, 20)*-2.0];
% y = [ones(1, 20)*-0.1 -0.1:-0.1:-2.0 ones(1, 20)*-2.0 -2.0:0.1:-0.1];

% x = [-4:0.1:-2 ones(1, 20)*-2 -2:-0.1:-4 ones(1, 20)*-4];
% y = [ones(1, 20)*-0.1 -0.1:-0.1:-2.1 ones(1, 20)*-2.1 -2.1:0.1:-0.1];

rx1 = [0, -2];
rx2 = [2, 0];
tx = [0, 0];
coor = [x(1), y(1)];
c = 1;
s = [];
t = linspace(0,2 * pi,1000);

for i = 1:length(x)
    if x(i) == 0
        continue
    end
    n_a1 = (sqrt((rx2(1)-x(i))^2 +(rx2(2)-y(i))^2)+sqrt((tx(1)-x(i))^2 +(tx(2)-y(i))^2)) / 2;
    n_b1 = sqrt(n_a1^2 - c^2);
    n_a2 = (sqrt((rx1(1)-x(i))^2 +(rx1(2)-y(i))^2)+sqrt((tx(1)-x(i))^2 +(tx(2)-y(i))^2)) / 2;
    n_b2 = sqrt(n_a2^2 - c^2);
    
%     e1_x = 1 + n_a1 * cos(t);
%     e1_y = n_b1 * sin(t);
%     e2_x = n_b2 * cos(t);
%     e2_y = -1 + n_a2 * sin(t);
%     plot(e1_x, e1_y);
%     hold on;
%     plot(e2_x, e2_y);
%     grid on;
%     disp(1);
%     pause(1);
%     close all;
    syms u v;
    eq1 = (u-1)^2/n_a1^2+v^2/n_b1^2-1;
    eq2 = u^2/n_b2^2+(v+1)^2/n_a2^2-1;
    s = [s solve(eq1,eq2,u,v)];
end

cach_u = [];
cach_v = [];
sel = [];

for i = 1:4
    if isreal(double(s(1).u(i)))
        sel = [sel i];
    end
end

for i = 1:length(s)
    for j = 1:4
        if isreal(double(s(i).u(j)))
            cach_u = [cach_u double(s(i).u(j))];
            cach_v = [cach_v double(s(i).v(j))];
        end
    end
end

scatter(cach_u, cach_v);
set(0,'defaultfigurecolor','w');
xlabel('x方向坐标/m');
ylabel('y方向坐标/m');
grid on;
