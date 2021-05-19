pathx = [ones(1, 10)*0.2 ones(1, 10)*0.3 ones(1, 10)*0.3 ones(1, 10)*0.2 ones(1, 10)*0.2];
pathy = [ones(1, 10)*(-0.4) ones(1, 10)*(-0.4) ones(1, 10)*(-0.3) ones(1, 10)*(-0.3) ones(1, 10)*(-0.4)];
% pathx = [ones(1, 10)*0.3 ones(1, 10)*0.2 ones(1, 10)*0.3 ones(1, 10)*0.3 ones(1, 10)*0.2];
% pathy = [ones(1, 10)*(-0.2) ones(1, 10)*(-0.3) ones(1, 10)*(-0.4) ones(1, 10)*(-0.2) ones(1, 10)*(-0.3)];
figure;
hold on;
dis = sqrt((pathx).^2 + (pathy).^2) + sqrt((pathx - 2 *ones(1, length(pathx))).^2 + (pathy).^2);
dis = dis / 0.333 * 2 * pi;
set(0,'defaultfigurecolor','w');
plot(dis, 'r','LineWidth',3);
title('position 1');
xlabel('time');
ylabel('phase');

figure;
hold on;
dis = sqrt((pathx).^2 + (pathy).^2) + sqrt((pathx).^2 + (pathy + 2 *ones(1, length(pathy))).^2);
dis = dis / 0.333 * 2 * pi;
set(0,'defaultfigurecolor','w');
plot(dis,'r','LineWidth',3);
title('position 2');
xlabel('time');
ylabel('phase');

