function [recover_signal,diff1, diff2] = CurveFitRecover(sin_signal,flip_idx,index_ls, chirp_samples, sample_rate, BW)
    recover_signal = sin_signal;
    for i = 1 : length(flip_idx)
        st = max(1,flip_idx(i) - 100);
        ed = min(flip_idx(i) + 100, length(recover_signal));
        recover_signal(st:ed) = zeros([1 (ed - st + 1)]);
    end
    diff1 = [];
    diff2 = [];

    last_diff = 1 + 0j;
    move_rate = 3e8 / 902e6;

    phase_change_threshold = 10 / move_rate * 2 * pi / sample_rate;
    for i = 1 : length(sin_signal) / chirp_samples

        st1 = chirp_samples * (i - 1) + 1;
        while(recover_signal(st1) == 0)
            st1 = st1 + 1;
        end
        ed1 = st1 + 1;
        while(recover_signal(ed1) ~= 0)
            ed1 = ed1 + 1;
        end
        ed1 = ed1 - 1;
        if(ed1 > chirp_samples * i)
            ed1 = chirp_samples * i;
        end
        
        ed2 = 0;
        st2 = ed1 + 1;
        while(recover_signal(st2) == 0)
            st2 = st2 + 1;
        end
        if st2 < chirp_samples * i
             ed2 = st2 + 1;
            while(recover_signal(ed2) ~= 0)
                ed2 = ed2 + 1;
            end
            ed2 = ed2 - 1;
            if(ed2 > chirp_samples * i)
                ed2 = chirp_samples * i;
            end
        end

        
        fun = @(x,xdata)x(1) + x(2)*xdata + x(3)*(xdata.^2);
        x0 = [1,1,1];
        if(ed1 - st1 > ed2 - st2)
            xdata = st1 - chirp_samples * (i - 1) : ed1 - chirp_samples * (i - 1);
            ydata = unwrap(angle(recover_signal(st1:ed1)));
            recover_signal(ed1 + 1:i * chirp_samples) = recover_signal(ed1 + 1:i * chirp_samples) * 0;
        else
            xdata = st2 - chirp_samples * (i - 1) : ed2 - chirp_samples * (i - 1);
            ydata = unwrap(angle(recover_signal(st2:ed2)));
            recover_signal((i - 1) * chirp_samples  + 1:st2 - 1) = recover_signal((i - 1) * chirp_samples  + 1:st2 - 1) * 0;
        end
        options = optimoptions('lsqcurvefit','Display', 'off');
        x = lsqcurvefit(fun,x0,xdata,ydata,[],[],options);
        rotate = ones([1 chirp_samples]);
        rotate(xdata) = exp(-1j* fun([0,x(2),x(3)],xdata));

        rotate = rotate ./ abs(rotate);
        recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) = recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) .* rotate;

%         if i > 9
%             figure;
%             set(0,'defaultfigurecolor','w');
%             plot(unwrap(angle(recover_signal(chirp_samples * (i - 8) + 1 : chirp_samples * i))));
%             xlabel('时间/采样点序列', 'FontSize',14);
%             ylabel('相位/弧度', 'FontSize',14);
%         end

        chirp_recover = recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i);
        zeros_idx_chirp = chirp_recover == 0;
        chirp_recover(zeros_idx_chirp) = [];
        chirp_mean = mean(chirp_recover);
        
        [cluster_idx,C_iq] = kmeans([real(chirp_recover)', imag(chirp_recover)'],2);      
        idx1 = cluster_idx == 1;
        idx2 = cluster_idx == 2;
        xdata1 = xdata(idx1);
        xdata2 = xdata(idx2);       
        ydata1 = ydata(idx1);
        ydata2 = ydata(idx2);
        x1 = lsqcurvefit(fun,x0,xdata1,ydata1,[],[],options);
        x2 = lsqcurvefit(fun,x0,xdata2,ydata2,[],[],options);
        rotate1 = exp(-1j* fun([0,x1(2),x1(3)],xdata1));
        iq1 = sin_signal(chirp_samples * (i - 1) + xdata1) .* rotate1;
        rotate2 = exp(-1j* fun([0,x2(2),x2(3)],xdata2));
        iq2 = sin_signal(chirp_samples * (i - 1) + xdata2) .* rotate2;
        
        C = [mean(iq1), mean(iq2)];
        chirp_diff = (C(1)-C(2)) * abs(chirp_mean) / chirp_mean;
        angle_test = imag(log(last_diff / chirp_diff));
        delta_idx = 10000000;
        if i > 1
            delta_idx = index_ls(i) - index_ls(i - 1);
        end
        diff1_phase = abs(imag(log(last_diff / ((C(1) - C(2)) * abs(C(2)) / C(2)))));
        diff2_phase = abs(imag(log(last_diff / ((C(2) - C(1)) * abs(C(1)) / C(1)))));
        if diff1_phase < diff2_phase && diff1_phase / (delta_idx) < phase_change_threshold
            diff1 = [diff1, (C(1) - C(2)) * abs(C(2)) / C(2)];
            diff2 = [diff2, (C(2) - C(1)) * abs(C(1)) / C(1)];

        elseif diff2_phase / (delta_idx) < phase_change_threshold
            diff1 = [diff1, (C(2) - C(1)) * abs(C(1)) / C(1)];
            diff2 = [diff2, (C(1) - C(2)) * abs(C(2)) / C(2)];
  
        else
            diff1 = [diff1, last_diff];
            diff2 = [diff2, -last_diff];
        end
        
        last_diff = diff1(end);
    end
     
end

