function [ydata] = FittedFunction(m,f_cfo,t_delta,flip, xdata, sample_rate, BW, flip_pos, chirp_samples)
    f0 = -BW / 2 + BW * (chirp_samples - flip_pos) / chirp_samples;
    k = BW / chirp_samples * sample_rate;
    ydata = zeros([1 length(xdata)]);
    for i = 1 : length(xdata)
        t = xdata(i) - 1;
        t = t / sample_rate;
        if abs(xdata(i) - flip_pos) <= 100
            ydata(i) = 0;
        elseif xdata(i) < flip_pos
            ydata(i) = 2 * pi * ((f0+f_cfo) * m * t_delta + 1/2*k*m*m*t_delta + (m*f_cfo+(m-1)*f0+k*m*m*t_delta) * t + 1/2*k*(m*m - 1)*t^2);
        else
            ydata(i) = 2 * pi * (flip + (f0-BW+f_cfo) * m * t_delta + 1/2*k*m*m*t_delta + (m*f_cfo+(m-1)*(f0-BW)+k*m*m*t_delta) * t + 1/2*k*(m*m - 1)*t^2);
        end        
    end
end