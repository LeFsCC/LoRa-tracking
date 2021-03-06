file_name = "square-y2";

diff_1 = load(strcat("../Data/5-22-result/", file_name, "-diff1.mat"));
diff_1 = diff_1.diff1;

diff_2 = load(strcat("../Data/5-22-result/", file_name, "-diff2.mat"));
diff_2 = diff_2.diff2;

diff_1 = unwrap(angle(diff_1));
diff_2 = unwrap(angle(diff_2));

window_width = 500;

figure;
set(0,'defaultfigurecolor','w');
plot(diff_1);
xlabel('采样序列/个', 'FontSize',14);
ylabel('相位/弧度', 'FontSize',14);

average_move = [];
step = 5;

for i = 1:step
   average_move = [average_move BreathFilter(diff_1((i-1)*length(diff_1)/step+1:i*length(diff_1)/step), window_width)];
end
% average_move = BreathFilter(diff_1, window_width);


figure;
set(0,'defaultfigurecolor','w');
plot(average_move,'r', 'LineWidth',1);
xlabel('采样序列/个', 'FontSize',14);
ylabel('相位/弧度', 'FontSize',14);
% average_move = average_move(1:100:length(average_move));
% intersection_tracking(average_move);

figure;
set(0,'defaultfigurecolor','w');
plot(diff_2);
% title('目标物体反射信号相位2');
xlabel('采样序列/个', 'FontSize',14);
ylabel('相位/弧度', 'FontSize',14);

f1 = strcat('../Data/5-22-result/', file_name, '-slide-diff1.mat');
save(f1, 'average_move');

average_move = [];

for i = 1:step
   average_move = [average_move BreathFilter(diff_2((i-1)*length(diff_2)/step+1:i*length(diff_2)/step), window_width)];
end
% average_move = BreathFilter(diff_2, window_width);

figure;
set(0,'defaultfigurecolor','w');
plot(average_move,'r', 'LineWidth',1);
xlabel('采样序列/个', 'FontSize',14);
ylabel('相位/弧度', 'FontSize',14);
% title('目标物体反射信号相位2');

f2 = strcat('../Data/5-22-result/', file_name, '-slide-diff2.mat');
save(f2, 'average_move');

% average_move = average_move(1:100:length(average_move));
% intersection_tracking(average_move);