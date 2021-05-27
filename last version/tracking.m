% pathx = [ones(1, 10)*0 ones(1, 10)*0.1 ones(1, 10)*0.1 ones(1, 10)*0 ones(1, 10)*0];
% pathy = [ones(1, 10)*(-0.1) ones(1, 10)*(-0.1) ones(1, 10)*(-0.2) ones(1, 10)*(-0.2) ones(1, 10)*(-0.1)];

pathx = [ones(1, 10)*0 ones(1, 10)*0.1 ones(1, 10)*0.1 ones(1, 10)*0];
pathy = [ones(1, 10)*(-0.2) ones(1, 10)*(-0.1) ones(1, 10)*(-0.2) ones(1, 10)*(-0.2)];

% pathx = [ones(1, 10)*0.3 ones(1, 10)*0.2 ones(1, 10)*0.3 ones(1, 10)*0.3 ones(1, 10)*0.2];
% pathy = [ones(1, 10)*(-0.2) ones(1, 10)*(-0.3) ones(1, 10)*(-0.4) ones(1, 10)*(-0.2) ones(1, 10)*(-0.3)];

figure;
hold on;
dis = sqrt((pathx).^2 + (pathy).^2) + sqrt((pathx - 2 *ones(1, length(pathx))).^2 + (pathy).^2);
dis = -dis / 0.332 * 2 * pi;
set(0,'defaultfigurecolor','w');
plot(dis, 'r','LineWidth',3);
title('横向接收机');
xlabel('时间');
ylabel('相位');

figure;
hold on;
dis = sqrt((pathx).^2 + (pathy).^2) + sqrt((pathx).^2 + (pathy + 2 *ones(1, length(pathy))).^2);
dis = -dis / 0.333 * 2 * pi;
set(0,'defaultfigurecolor','w');
plot(dis,'r','LineWidth',3);
title('纵向接收机');
xlabel('时间');
ylabel('相位');


figure;
scatter(pathx, pathy);

