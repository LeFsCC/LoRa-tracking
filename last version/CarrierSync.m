function [signal_sync] = CarrierSync(signal_input, BW, SF, sample_rate, pkt_size)
    signal_sync = signal_input;
    chirp_samples = 2^SF * sample_rate / BW ;
    pkt_number = floor(length(signal_input) / pkt_size / chirp_samples);

    for i = 1 : pkt_number
        pkt_start = 1 + (i - 1) * pkt_size * chirp_samples;
        pkt_end = i * pkt_size * chirp_samples;
        signal = signal_sync(pkt_start: pkt_end);
        angles = unwrap(angle(signal));
        slide_mean = zeros(size(angles));
        window_size = 50;
        for j = 1: length(angles) - window_size
            slide_mean(j) = mean(angles(j:j + window_size));
        end
        delta_angle = (slide_mean(chirp_samples - window_size) - slide_mean(5)) / (chirp_samples - window_size - 5);
        signal_sync(pkt_start: pkt_end) = signal_sync(pkt_start: pkt_end) .* exp([1:length(signal)] * delta_angle * (-1j)); 
        %% downchirp ÖØÐÂ²¹³¥
        signal = signal_sync(pkt_start: pkt_end);
        angles = unwrap(angle(signal(1 + 10 * chirp_samples:chirp_samples * 2 + 10 * chirp_samples)));
        slide_mean = zeros(size(angles));
        window_size = 50;
        for j = 1: length(angles) - window_size
            slide_mean(j) = mean(angles(j:j + window_size));
        end
        delta_angle = (slide_mean(chirp_samples - window_size) - slide_mean(1)) / (chirp_samples - window_size - 1);
        signal_sync(pkt_start + 10 * chirp_samples + 1 : pkt_start + 12.25 * chirp_samples) = signal_sync(pkt_start + 10 * chirp_samples + 1 : pkt_start + 12.25 * chirp_samples) .* exp([1:2.25 * chirp_samples] * delta_angle * (-1j)); 

    end

%     angles = unwrap(angle(signal_input(1:chirp_samples * 2 - 1)));
%     slide_mean = zeros(size(angles));
%     window_size = 50;
%     for i = 1: length(angles) - window_size
%         slide_mean(i) = mean(angles(i:i + window_size));
%       
%     end
%     delta_angle = (slide_mean(chirp_samples - window_size) - slide_mean(1)) / (chirp_samples - window_size - 1);
%     signal_sync = signal_sync .* exp([1:length(signal_sync)] * delta_angle * (-1j));
%     
%     angles = unwrap(angle(signal_sync));
% 
%         
%     figure;
%     plot(angles(1:chirp_samples * 20));
end