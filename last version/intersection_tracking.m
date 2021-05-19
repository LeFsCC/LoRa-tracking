function [] = intersection_tracking(average_move)
    rx = [0, 0];
    tx = [3, 0];
    coor = [tx(1)-1.5, -0.4];
    a = (sqrt((rx(1)-coor(1))^2 +(rx(2)-coor(2))^2)+sqrt((tx(1)-coor(1))^2 +(tx(2)-coor(2))^2)) / 2;
    c = (rx(1) - tx(1)) / 2;
    b = sqrt(a^2 - c^2);
    wave_length = 3e8 / 902e6;
    phase_shift = average_move(1) - 0;
    s = [];
    ind = 1:1:length(average_move)-1;
    for i = 2:length(average_move)
        phase_shift = average_move(i) - phase_shift;
        a = a + phase_shift * wave_length / (2 * pi);
        b = sqrt(a^2 - c^2);
        syms x;
        eq=(x-1.5)^2/a^2+0.4^2/b^2-1;
        s=[s solve(eq,x)];
        phase_shift = average_move(i);
    end
    figure;
    plot(ind, s(1, :));
    title('目标物体运动追踪');
    xlabel('时间序列');
    ylabel('距离');
    figure;
    plot(ind, s(2, :));
    title('目标物体运动追踪');
    xlabel('chirp number');
    ylabel('phase/rad');
end

