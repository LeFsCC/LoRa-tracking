function [average_move] = BreathFilter(diff)
    move_data = unwrap(angle(diff));
    slide_window = 500;
    move_interp = move_data;
    average_move = zeros([1 (length(move_interp) - slide_window)]);
    for i = 1:length(average_move)
        average_move(i) = mean(move_interp(i:i + slide_window));
    end
    for i = 1:length(average_move)
        average_move(i) = mean(move_interp(i:i + slide_window));
    end
    figure;
    set(0,'defaultfigurecolor','w');
    plot(average_move,'r', 'LineWidth',1);
    xlabel('chirp number');
    ylabel('phase');
end