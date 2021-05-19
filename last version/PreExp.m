BW = 500000;
SF = 11;
sample_rate = 1e6;
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
pkt_size = 12.25 + 18;
pkt_interval = 0.3;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');

%fid = fopen("Data/commercial/1-4/1_4_zoulang_wuren",'rb');
fid = fopen("Data/210301/huxi_3",'rb');

raw_data =fread(fid,[2 inf],'float');

fclose(fid);
raw_data_complex = complex(raw_data(1,:),raw_data(2,:));
figure;
plot(abs(raw_data_complex));

% figure;
% pspectrum(raw_data_complex(280000:450000), sample_rate, 'spectrogram');
[data, flip_idx, position, sin_signal, time_ls] = PacketContentDetect(raw_data_complex,BW,SF,sample_rate,pkt_size, pkt_interval);
figure;
plot(unwrap(angle(raw_data_complex(position(2):position(2) + chirp_samples * 30))));
figure;
plot(unwrap(angle(data(position(2):position(2) + chirp_samples * 30))));

figure;
hold on;
temp_data = raw_data_complex(position(2):position(2) + chirp_samples - 1);
temp_mean = ones([1 chirp_samples]) * mean(temp_data);

plot(real(temp_data),'Color','#1F9137');
plot(real(temp_mean),'-*','Color','#69E082');
plot(imag(temp_data),'Color','#942835');
plot(imag(temp_mean),'-*','Color','#FF8593');


% featureChirp =[repmat(up_chirp,1,1),zeros([1, 2 * chirp_samples]),repmat(down_chirp,1,2), down_chirp(1:chirp_samples/4)];
% ploting = data(position(2): position(2) + pkt_size * 2 * chirp_samples);
% [corr_result, lag] = xcorr(ploting ./ abs(ploting), featureChirp);
% corr_result = corr_result / length(featureChirp);
% figure;
% plot(abs(real(corr_result)));
% figure;
% hold on;
% plot(real(data(position(2): position(2) + pkt_size * 2 * chirp_samples)));
% plot(imag(data(position(2): position(2) + pkt_size * 2 * chirp_samples)));
figure;
plot(abs(data(position(2): position(2) + pkt_size * 2 * chirp_samples)));

for i = 1 : length(position)
    if mod(i, 25) == 1
        figure;
    end
    plot_idx = mod(i,25);
    if plot_idx == 0
        plot_idx = 25;
    end
    subplot(5,5,plot_idx);
    plot(unwrap(angle(data(position(i):position(i) + chirp_samples * 30))));
    %plot(abs(data(position(i):position(i) + chirp_samples * 30)));
end

figure;
plot(unwrap(angle(data(position(9):position(9) + chirp_samples * 30))));
figure;
plot(abs(data(position(9):position(9) + chirp_samples * 30)));

% figure;
% plot(unwrap(angle(data(position(17):position(17) + pkt_size * chirp_samples))));
% signal_chirp_size = pkt_size - 1.25;
% figure;
% for i = 2 : length(position)
%     st = chirp_samples * signal_chirp_size * (i - 1) + chirp_samples + 1;
%     lg =    1 * chirp_samples;
%     %lg = length(sin_signal) - 1;
%     signal_scatter = sin_signal(st:1:st + lg);
%     % signal_scatter = sin_signal(st:1:end);
%     c = linspace(1,10,length(signal_scatter));
%     subplot(5,4,i);
% %     scatter3(real(signal_scatter), imag(signal_scatter), 1:length(signal_scatter), 10,c);
%     scatter(real(signal_scatter), imag(signal_scatter), 10,c);
%     xlabel('real');
%     ylabel('imag');
%     zlabel('time');
% end
% 
% pkt_idx = 3;
% figure;
% for i = 1:signal_chirp_size
%     st = chirp_samples * signal_chirp_size * (pkt_idx-1) + (i - 1)  * chirp_samples + 1;
%     lg = 1 * chirp_samples;
%     signal_scatter = sin_signal(st:1:st + lg);
%     c = linspace(1,10,length(signal_scatter));
%     subplot(6,6,i);
%     scatter(real(signal_scatter), imag(signal_scatter), 10,c);
%     xlabel('real');
%     ylabel('imag');
%     zlabel('time');
% end
% 
% 
% chirp_idx = 2;
% figure;
% st = chirp_samples * signal_chirp_size * (pkt_idx-1) + (chirp_idx - 1)  * chirp_samples + 1;
% lg = 17 * chirp_samples;
% signal_scatter = sin_signal(st:1:st + lg);
% c = linspace(1,10,length(signal_scatter));
% scatter3(real(signal_scatter), imag(signal_scatter),1:length(signal_scatter), 10,c);
% xlabel('real');
% ylabel('imag');
% zlabel('time');
% 
% figure;
% plot(unwrap(angle(signal_scatter)));
% figure;
% plot(abs((signal_scatter)));


% sfr = StraightFlipRecover(sin_signal,flip_idx,chirp_samples);
% sfr_diff = sfr(2:end) - sfr(1: end - 1);
% diff_scatter = sfr_diff;
% figure;
% c = linspace(1,10,length(diff_scatter));
% scatter3(real(diff_scatter), imag(diff_scatter),1:length(diff_scatter), 1,c);
% xlabel('real');
% ylabel('imag');
% zlabel('time');
% 
% figure;
% plot(abs(sfr_diff));
[cfc, diff1, diff2] = CurveFitRecover(sin_signal,flip_idx,chirp_samples, sample_rate, BW);
normal_diff1 = diff1 ./ abs(diff1);
c = linspace(1,10,length(normal_diff1));
figure;
scatter3(real(diff1), imag(diff1), 1:length(diff1),20,c);
figure;
scatter3(real(normal_diff1), imag(normal_diff1), 1:length(normal_diff1),20,c);
figure;
plot(unwrap(angle(diff1)));

% normal_diff2 = diff2 ./ abs(diff2);
% c = linspace(1,10,length(normal_diff2));
% figure;
% scatter3(real(diff2), imag(diff2), 1:length(diff2),20,c);
% figure;
% scatter3(real(normal_diff2), imag(normal_diff2), 1:length(normal_diff2),20,c);
% figure;
% plot(unwrap(angle(diff2)));
% diff_envelope = [0];
% % for i = 1 : length(diff) / 50 - 1
% %     [M1,idx1] = max(abs(diff(1 + (i - 1) * 50 : 1 + (i - 1) * 50 +25)));
% %     [M2,idx2] = max(abs(diff(1 + (i - 1) * 50 + 25 : 50 * i)));
% %     flip1 = diff((i - 1) * 50 + idx1);
% %     flip2 = diff((i - 1) * 50 + 25 + idx2);
% %     angle_test1 = imag(log(diff_envelope(end) / flip1));
% %     angle_test2 = imag(log(diff_envelope(end) / flip2));
% %     if(abs(angle_test1) < abs(angle_test2))
% %         diff_envelope = [diff_envelope, diff((i - 1) * 50 + idx1)];
% %     else
% %         diff_envelope = [diff_envelope, diff((i - 1) * 50 + 25 + idx2)];
% %     end
% % end
% for i = 1 : length(diff) / 25 - 1
%     [M,idx] = max(abs(diff(1 + (i - 1) * 25 : 1 + i * 25)));
%     flip = diff((i - 1) * 25 + idx);
%     angle_test = imag(log(diff_envelope(end) / flip));
%     if(abs(angle_test) < 0.5*pi)
%         diff_envelope = [diff_envelope, diff((i - 1) * 25 + idx)];
%     end
% end
% % diff_envelope = diff_envelope ./ abs(diff_envelope);
% diff_scatter = diff_envelope(1 : end);
% figure;
% c = linspace(1,10,length(diff_scatter));
% scatter3(real(diff_scatter), imag(diff_scatter),1:length(diff_scatter), 10,c);
% xlabel('real');
% ylabel('imag');
% zlabel('time');
% 
% figure;
% plot(unwrap(angle(diff_scatter)));
% 
% figure;
% plot(abs(diff));


x1 = 20;
x2 = 80;
y = 1;
delta = 0.1;
d = sqrt(x1*x1+y*y) + sqrt(x2*x2 + y * y)
sqrt((x1-delta)*(x1-delta)+y*y) + sqrt((x2+delta)*(x2+delta) + y * y) - d
sqrt((x1+delta)*(x1+delta)+y*y) + sqrt((x2-delta)*(x2-delta) + y * y) - d
sqrt(x1*x1+(y-delta)*(y-delta)) + sqrt(x2*x2 + (y-delta) * (y-delta)) - d
sqrt(x1*x1+(y+delta)*(y+delta)) + sqrt(x2*x2 + (y+delta) * (y+delta)) - d
