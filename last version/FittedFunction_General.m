function [ydata] = FittedFunction_General(x, xdata, flip_pos)
    ydata = zeros([1 length(xdata)]);
    for i = 1 : length(xdata)
        if abs(xdata(i) - flip_pos) <= 100
            ydata(i) = 0;
        elseif xdata(i) < flip_pos
            ydata(i) = x(1) + x(2) * xdata(i) + x(3)*xdata(i)^2;
        else
            ydata(i) = x(1) + x(4) + (x(2)+x(5)) * xdata(i) + x(3)*xdata(i)^2;
        end        
    end
end