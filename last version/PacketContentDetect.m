function [data,idx_ls, position, sin_signal, index_ls] = PacketContentDetect(raw_data_complex, BW, SF, sample_rate, pkt_size)

chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');
featureChirp =[repmat(up_chirp,1,2),repmat(down_chirp,1,2), down_chirp(1:chirp_samples/4)];

%% 手动修正CFO

% [corr_result, lag] = xcorr(raw_data_complex, featureChirp);
% corr_result = corr_result / length(featureChirp);
% [~,I] = max(abs(corr_result));
% CFO_detect_idx = lag(I) - 2 * chirp_samples;
% data_phase = unwrap(angle(raw_data_complex));
% CFO_rate = (data_phase(CFO_detect_idx + 2 * chirp_samples) - data_phase(CFO_detect_idx)) / 2 / chirp_samples;

CFO_rate = 0.074992954799107;

% CFO_rate = 0.077;

data = raw_data_complex .*  exp(-1j * (1:length(raw_data_complex)) * CFO_rate);

featureChirp =[repmat(up_chirp,1,2),repmat(down_chirp,1,2), down_chirp(1:chirp_samples/4)];

[corr_result, lag] = xcorr(data, featureChirp);

corr_result = corr_result / length(featureChirp);

position = [];

corr_threshold = 0.05;
abs_threshold = 0.1;
i = length(data) + 1 + 12.25 * chirp_samples - length(featureChirp);
pkt_start_compensate = 8 * chirp_samples;
guard = 300;

while (i <= length(corr_result) -  pkt_size * chirp_samples)
    [M,I] = max(abs(corr_result(i: i + guard * 2)));
    if I == 1
        if abs(corr_result(i)) < corr_threshold  || abs(data(lag(i) - 1 - pkt_start_compensate + guard)) < abs_threshold || abs(data(lag(i) - 1 - pkt_start_compensate + pkt_size * chirp_samples - guard)) < abs_threshold
            i = i + 1;
            continue;
        end
        if min(abs(data(lag(i) - 1 - pkt_start_compensate + guard:lag(i) - 1 - pkt_start_compensate + pkt_size * chirp_samples - guard))) < abs_threshold / 2
            i = i + 1;
            continue;
        end
        position = [position (lag(i) - pkt_start_compensate)];
        i = i + pkt_size * chirp_samples - 1;
    else
        i = i + I - 2; 
    end
    i = i + 1;
end

%% 修正xcorr offset
sin_signal = [];
idx_ls  = [1];
signal_chirp_size = pkt_size - 1.25;
position = position(2:end);
index_ls = [];

for i = 1 : length(position)
    y = fft(data(position(i) + chirp_samples: position(i) + chirp_samples * 2 - 1).* down_chirp);
    [M,idx] = max(abs(y));
    if idx > 2048
        idx = idx - 4096;
    end
    offset = [2 * (1 - idx) - 1, 2 * (1 - idx), 2 * (1 - idx) + 1];
    m_candidate = [0,0,0];
    for j = 1:3
        temp = abs(fft(data(position(i) + chirp_samples + offset(j): position(i) + chirp_samples * 2 - 1 + offset(j)).* down_chirp)); 
        m_candidate(j) = temp(1);
    end
    [~, idx] = max(m_candidate);
%     position(i) = position(i) + offset(idx);
    %% fft decode
    pkt_start = (i - 1) * signal_chirp_size * chirp_samples + 1;
    pkt_conj = [];
    start_chirp = 2;
    for j = start_chirp : 10
        idx_ls = [idx_ls (pkt_start + (j - 1) * chirp_samples)];

        y = fft(data(position(i) + (j - 1) * chirp_samples: position(i) + j * chirp_samples - 1).* down_chirp);
        [M,idx] = max(abs(y));

        if idx > 2048
            idx = idx - 2048;
        end
        f_start = BW / 2 - (idx - 1) * sample_rate / chirp_samples;
        f_k = -BW / chirp_duration;
        t_flip = (- BW / 2 - f_start) / f_k;
        conj_chirp = [chirp(1/sample_rate:1/sample_rate:t_flip, f_start, chirp_duration, f_start - BW, 'linear',0,'complex') chirp(t_flip+ 1/sample_rate:1/sample_rate:chirp_duration, f_start + BW, chirp_duration, f_start, 'linear',0,'complex')];
        pkt_conj = [pkt_conj conj_chirp];
        if(idx ~= 1)
            idx_ls = [idx_ls (pkt_start + (j - 2) * chirp_samples + (4098 - idx * 2))];
        end
    end
    pkt_conj = [pkt_conj, up_chirp, up_chirp];
    sin_signal = [sin_signal data(position(i) + chirp_samples:position(i) + 12 * chirp_samples - 1).*pkt_conj];
    index_ls = [index_ls ((position(i)+chirp_samples):chirp_samples:(position(i) + 12 * chirp_samples - 1))];
    idx_ls = [idx_ls (pkt_start + (11 - 1) * chirp_samples) (pkt_start + (12 - 1) * chirp_samples)];
    
%     figure;
%     plot(unwrap(angle(data(position(i) + chirp_samples:position(i) + 12 * chirp_samples - 1).*pkt_conj)));
%     set(0,'defaultfigurecolor','w');
%     xlabel('时间/采样点序列', 'FontSize',14);
%     ylabel('相位/弧度', 'FontSize',14);
%     
%     dechirp = data(position(i) + chirp_samples:position(i) + 12 * chirp_samples - 1).*pkt_conj;
    
%     figure;
%     Y = abs(fft(dechirp(1:chirp_samples)));
%     plot(Y())
%     set(0,'defaultfigurecolor','w');
%     xlabel('FFT/bin', 'FontSize',14);
%     ylabel('频率/Hz', 'FontSize',14);
    
    pkt_conj = [];
    
    for j = 1 : pkt_size - 12.25
        idx_ls = [idx_ls (pkt_start + (12 + j - 1) * chirp_samples)]; % 此处加的点是chirp尾

        y = fft(data(position(i) + (j - 1 + 12.25) * chirp_samples: position(i) + (j + 12.25) * chirp_samples - 1).* down_chirp);
        [M,idx] = max(abs(y));
        if idx > 2048
            idx = idx - 2048;
        end
        if(idx ~= 1)
            idx_ls = [idx_ls (pkt_start + (12 + j - 2) * chirp_samples + (4098 - idx * 2))];
        end
        f_start = BW / 2 - (idx - 1) * sample_rate / chirp_samples;
        f_k = -BW / chirp_duration;
        t_flip = (- BW / 2 - f_start) / f_k;
        conj_chirp = [chirp(0:1/sample_rate:t_flip - 1/sample_rate, f_start, chirp_duration, f_start - BW, 'linear',0,'complex') chirp(t_flip:1/sample_rate:chirp_duration - 1/sample_rate, f_start + BW, chirp_duration, f_start, 'linear',0,'complex')];
        pkt_conj = [pkt_conj conj_chirp];
    end
    
    sin_signal = [sin_signal data(position(i) + 12.25 * chirp_samples:position(i) + pkt_size * chirp_samples - 1).*pkt_conj];
    index_ls = [index_ls ((position(i) + 12.25 * chirp_samples):chirp_samples:(position(i) + pkt_size * chirp_samples - 1))];
end