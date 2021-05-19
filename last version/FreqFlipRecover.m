% 频率跳变的复原，基于2分类
function [signal_output] = FreqFlipRecover(signal_input, sample_per_chirp, pkt_size, flip_offset)
    phase_input = unwrap(angle(signal_input));
    signal_output = signal_input;
    pkt_number = floor(length(phase_input) / pkt_size / sample_per_chirp);
    last_level = 1;
    last_sample = 1 + 1j * 0;
    for i = 1: pkt_number % 每个packet
        pkt_start = 1 + (i - 1) * pkt_size * sample_per_chirp;
        for j = 1: 12 % 12个完整的up/down chirp，跳变发生在chirp结束时
            chirp_start = pkt_start + (j - 1) * sample_per_chirp;
            chirp_end = chirp_start - 1 + sample_per_chirp;
            [signal_output(chirp_start:chirp_end), last_sample, last_level] = CheckDif(signal_input(chirp_start:chirp_end), last_sample, last_level);
        end
        %0.25个down chirp，跳变发生在结束时
        [signal_output(pkt_start + 12 * sample_per_chirp: pkt_start - 1 + 12.25 * sample_per_chirp), last_sample, last_level] = CheckDif(signal_input(pkt_start + 12 * sample_per_chirp: pkt_start - 1 + 12.25 * sample_per_chirp), last_sample, last_level);
        for j = 1 : length(flip_offset) % data chirp， 跳变位置发生在chirp之间，以及flip_offset(j)处
            chirp_start = pkt_start + (j - 1 + 12.25) * sample_per_chirp;
            chirp_flip = chirp_start - 1 + flip_offset(j);
            chirp_end = chirp_start - 1 + sample_per_chirp;
            [signal_output(chirp_start:chirp_flip), last_sample, last_level] = CheckDif(signal_input(chirp_start:chirp_flip), last_sample, last_level);
            [signal_output(chirp_flip + 1:chirp_end), last_sample, last_level] = CheckDif(signal_input(chirp_flip + 1:chirp_end), last_sample, last_level);
        end
    end
end

function [signal_output, last_sample, last_level] = CheckDif(signal_input,last_sample, last_level)
    idx = kmeans([real(signal_input)', imag(signal_input)'],2);
    check_idx = 1;
    count = 0;
    while count < 3
        if idx(check_idx) == last_level
            count = count + 1;
        else
            count = 0;
        end
        check_idx = check_idx + 1;
    end
    ref = mean(signal_input(check_idx - 3:check_idx - 1));
    dif = last_sample / ref / abs(last_sample / ref);
    signal_output = signal_input * dif;
    last_level = idx(end);
    last_sample = signal_output(end);
end