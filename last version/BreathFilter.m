function [average_move] = BreathFilter(move_data, window_width)
    slide_window = window_width;
    move_interp = move_data;
    average_move = zeros([1 (length(move_interp) - slide_window)]);
    for i = 1:length(average_move)
        average_move(i) = mean(move_interp(i:i + slide_window));
    end
    for i = 1:length(average_move)
        average_move(i) = mean(move_interp(i:i + slide_window));
    end
end