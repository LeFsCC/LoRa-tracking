
down_chirp = chirp(0:0.01:2-0.01, 0, 1.99, 15,'linear',0);
set(0,'defaultfigurecolor','w')
plot(down_chirp);
xlabel('时间/秒');
ylabel('幅度');

