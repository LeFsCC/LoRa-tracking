function [signal_output] = OutlierRecover(signal_input, sample_per_chirp, pkt_size)
    phase_input = unwrap(angle(signal_input));
    signal_output = signal_input;
    delta_phase = phase_input - [phase_input(1) phase_input(1:end - 1)];
    %outlier_rate = [0,0, abs(delta_phase(3:end) ./ delta_phase(2:end - 1)); abs(delta_phase(2:end - 1) ./ delta_phase(3:end)), 0, 0];
% figure;
% hold on;
% plot(outlier_rate(1,:));
% plot(outlier_rate(2,:));
    delta_mean = zeros(1, length(delta_phase));
    pkt_number = floor(length(delta_phase) / pkt_size / sample_per_chirp);
    for i = 1: pkt_number
        pkt_start = 1 + (i - 1) * pkt_size * sample_per_chirp;
        for j = 1: 12
            chirp_start = pkt_start + (j - 1) * sample_per_chirp;
            chirp_end = chirp_start - 1 + sample_per_chirp;
            delta_mean(chirp_start: chirp_end) = ones(1, sample_per_chirp) * mean(abs(delta_phase(chirp_start + 3 : chirp_end - 3)));
        end
        delta_mean(pkt_start + 12 * sample_per_chirp: pkt_start - 1 + 12.25 * sample_per_chirp) = ones(1, sample_per_chirp * 0.25) * mean(abs(delta_phase(pkt_start + 12 * sample_per_chirp + 3: pkt_start - 1 + 12.25 * sample_per_chirp - 3)));
        for j = 1 : (pkt_size - 12.25)
            chirp_start = pkt_start + (j - 1 + 12.25) * sample_per_chirp;
            chirp_end = chirp_start - 1 + sample_per_chirp;
            delta_mean(chirp_start: chirp_end) = ones(1, sample_per_chirp) * mean(abs(delta_phase(chirp_start + 3 : chirp_end - 3)));
        end
    end
    outlier_rate = [0, abs(delta_phase(2:end) ./ delta_mean(2:end));abs(delta_phase(2:end) ./ delta_mean(2:end)), 0];
    
    outlier_threshold = 2;
    for i = 1:length(signal_output)
        if(outlier_rate(1,i) > outlier_threshold && outlier_rate(2,i) > outlier_threshold)
            signal_output(i) = signal_output(i - 1);
            continue;
        end
    end
end