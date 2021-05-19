function [signal_output] = FlipRecover(signal_input, sample_per_chirp, pkt_size, flip_offset)
    phase_input = unwrap(angle(signal_input));
    signal_output = signal_input;
    phase_mean = zeros(1, length(phase_input));
    pkt_number = floor(length(phase_input) / pkt_size / sample_per_chirp);
    for i = 1: pkt_number
        pkt_start = 1 + (i - 1) * pkt_size * sample_per_chirp;
        for j = 1: 12
            chirp_start = pkt_start + (j - 1) * sample_per_chirp;
            chirp_end = chirp_start - 1 + sample_per_chirp;
            phase_mean(chirp_start: chirp_end) = ones(1, sample_per_chirp) * mean(phase_input(chirp_start + 3 : chirp_end - 3));
        end
        phase_mean(pkt_start + 12 * sample_per_chirp: pkt_start - 1 + 12.25 * sample_per_chirp) = ones(1, sample_per_chirp * 0.25) * mean(phase_input(pkt_start + 12 * sample_per_chirp + 3: pkt_start - 1 + 12.25 * sample_per_chirp - 3));
        for j = 1 : length(flip_offset)
            chirp_start = pkt_start + (j - 1 + 12.25) * sample_per_chirp;
            chirp_flip = chirp_start - 1 + flip_offset(j);
            chirp_end = chirp_start - 1 + sample_per_chirp;
            
            phase_mean(chirp_start: chirp_flip) =  ones(1, flip_offset(j)) * mean(phase_input(chirp_start + 3 : chirp_flip - 3));
            phase_mean(chirp_flip + 1: chirp_end) =  ones(1, sample_per_chirp - flip_offset(j)) * mean(phase_input(chirp_flip + 3 : chirp_end - 3));
        end
    end
%     chirp_number =  floor(length(phase_input) / sample_per_chirp);
% 
%     for i = 1: chirp_number
%         phase_mean((i - 1) * sample_per_chirp + 1 : i * sample_per_chirp) = ones(1, sample_per_chirp) * mean(phase_input((i - 1) * sample_per_chirp + 5 : i * sample_per_chirp-5));
%     end
%     figure;
%     plot(phase_mean);
    for i = 1:length(signal_output)
        signal_output(i) = signal_input(i) .* exp(1j * (- phase_mean(i)));
    end
end