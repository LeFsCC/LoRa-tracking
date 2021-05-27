file_name = "square-y3";

average_move_y1 = load(strcat("../Data/5-22-result/", file_name, "-slide-diff2.mat"));
average_move_y1 = average_move_y1.average_move;

file_name = "square-x8";
average_move_x = load(strcat("../Data/5-22-result/", file_name, "-slide-diff1.mat"));
average_move_x = average_move_x.average_move;

step = 5;

average_move_y = [];
for i=1:step
    average_move_y = [average_move_y average_move_y1((i-1)*length(average_move_y1)/step+400:(i-1)*length(average_move_y1)/step+length(average_move_x)/step+399)];
end

tx = [0, 0];
rx1 = [0, -2]; % 对应y
rx2 = [2, 0];  % 对应x
st_pos = [0, -0.1];
c = 1;

phase_diff_y = 0;
phase_diff_x = 0;

n_a1 = (sqrt((rx1(1)-st_pos(1))^2 +(rx1(2)-st_pos(2))^2)+sqrt((tx(1)-st_pos(1))^2 +(tx(2)-st_pos(2))^2)) / 2;
n_b1 = sqrt(n_a1^2 - c^2);
n_a2 = (sqrt((rx2(1)-st_pos(1))^2 +(rx2(2)-st_pos(2))^2)+sqrt((tx(1)-st_pos(1))^2 +(tx(2)-st_pos(2))^2)) / 2;
n_b2 = sqrt(n_a2^2 - c^2);
s = [];

average_move_x = average_move_x(1:10:end);
average_move_y = average_move_y(1:10:end);

figure;
plot(average_move_x);

figure;
plot(average_move_y);

for i = 2:length(average_move_x)
    phase_diff_y = (average_move_y(i) - average_move_y(1)) * 0.333 / (2 * pi);
    phase_diff_x = (average_move_x(i) - average_move_x(1)) * 0.333 / (2 * pi);
    a1 = n_a1 - phase_diff_y;
    b1 = sqrt(a1^2 - c^2);
    a2 = n_a2 - phase_diff_x;
    b2 = sqrt(a2^2 - c^2);

    syms u v;
    eq1 = (u-c)^2/a2^2+v^2/b2^2-1;
    eq2 = u^2/b1^2+(v+c)^2/a1^2-1;
    s = [s solve(eq1,eq2,u,v)];
end

cach_u = [];
cach_v = [];

for i = 1:length(s)
    for j = 1:4
        if isreal(double(s(i).u(j)))
            cach_u = [cach_u double(s(i).u(j))];
            cach_v = [cach_v double(s(i).v(j))];
        end
    end
end

cach_u1 = cach_u(1:2:end);
cach_v1 = cach_v(1:2:end);

cach_u2 = cach_u(2:2:end);
cach_v2 = cach_v(2:2:end);

figure;
set(0,'defaultfigurecolor','w');
scatter(cach_u1, cach_v1);
xlabel('横坐标/m', 'FontSize',14);
ylabel('纵坐标/m', 'FontSize',14);
hold on;
scatter(cach_u2, cach_v2);
xlabel('横坐标/m', 'FontSize',14);
ylabel('纵坐标/m', 'FontSize',14);

% temp_u = [];
% temp_v = [];
% 
% for i = 1:length(cach_v)
%     if cach_v(i) < 0
%         temp_u = [temp_u cach_u(i)];
%         temp_v = [temp_v cach_v(i)];
%     end
% end
% 
% figure;
% scatter(temp_u, temp_v);
% xlabel('横坐标/m', 'FontSize',14);
% ylabel('纵坐标/m', 'FontSize',14);
% set(0,'defaultfigurecolor','w');





