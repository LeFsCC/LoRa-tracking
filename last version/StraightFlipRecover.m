function [recover_signal] = StraightFlipRecover(sin_signal,flip_idx, chirp_samples)
    recover_signal = sin_signal;
    for i = 1 : length(flip_idx)
        st = max(1,flip_idx(i) - 100);
        ed = min(flip_idx(i) + 100, length(recover_signal));
        recover_signal(st:ed) = zeros([1 (ed - st + 1)]);
    end
    x=  recover_signal == 0;
    figure;
    scatter3(real(recover_signal(100:400)), imag(recover_signal(100:400)), 1:length(recover_signal(100:400)), 10);
    last = 1 + 0j;
    for i = 1 : length(sin_signal) / 4096
        
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
        rotate = recover_signal(ed1) / recover_signal(st1);
        rotate = rotate / abs(rotate);
        count = ed1 - st1;
        
       
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
            count = count + ed2 - st2;
            rotate2 = recover_signal(ed2) / recover_signal(st2);
            rotate2 = rotate2 / abs(rotate2);
            rotate = rotate * rotate2;
            
        end
        rotate_rate = rotate^(1/count);
        recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) = recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) .* ( (rotate_rate).^((0:chirp_samples - 1)));
        flip_rotate = recover_signal(st1) / last;
        flip_rotate = flip_rotate / abs(flip_rotate);
        recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) = recover_signal(chirp_samples * (i - 1) + 1 : chirp_samples * i) / flip_rotate;
        if(ed2 > 0)
            flip_rotate = recover_signal(st2) / recover_signal(ed1);
            flip_rotate = flip_rotate / abs(flip_rotate);
            recover_signal(st2 : chirp_samples * i) = recover_signal(st2 : chirp_samples * i) / flip_rotate;
            last = recover_signal(ed2);
        else
            last = recover_signal(ed1);
        end
    end
    recover_signal(x) = [];
%     figure;
%     plot(unwrap(angle(recover_signal)));
    
    figure;
    scatter3(real(recover_signal(1:chirp_samples * 1)), imag(recover_signal(1:chirp_samples * 1)), 1:length(recover_signal(1:chirp_samples * 1)), 1);
    figure;
    plot(unwrap(angle(recover_signal(1:chirp_samples * 1))));
end