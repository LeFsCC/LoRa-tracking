
BW = 125000;
SF = 7;
sample_rate = 1e6;
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
pkt_size = 29.25;

fid = fopen("Data/1-14/1-14",'rb');
% fid = fopen("Data/lora_sf7_payload64",'rb');
%fid = fopen("Data/lora_sf7",'rb');
raw_data =fread(fid,[2 inf],'float');
fclose(fid);
data_tx = complex(raw_data(1,:),raw_data(2,:));
figure;
hold on;
plot(real(data_tx));
plot(imag(data_tx));
figure;
hold on;
st = 0.7e6;
temp_data = data_tx(st:st+chirp_samples - 1);
temp_mean = ones([1 chirp_samples]) * mean(temp_data);

plot(real(temp_data),'Color','#1F9137');
plot(real(temp_mean),'-*','Color','#69E082');
plot(imag(temp_data),'Color','#942835');
plot(imag(temp_mean),'-*','Color','#FF8593');


simu_data = repmat(data_tx, 1, 10);
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');
figure;
plot(unwrap(angle(data_tx)));

temp_inspector = 12.25 * 1024 + 1;
figure;
hold on;
plot(real(data_tx(temp_inspector:temp_inspector+1023)));
plot(imag(data_tx(temp_inspector:temp_inspector+1023)));


% temp_inspector = 10473 + 1024*0 + 4;
temp_inspector = position(2) + 12.25 * chirp_samples;
temp_inspector_end = temp_inspector + 2 * 1024  * 1;
% data_norm = data ./ abs(data);
temp_data = data(temp_inspector:temp_inspector_end);
% figure;
% pspectrum(temp_data, sample_rate, 'spectrogram');
figure;
hold on;
plot(real(data(temp_inspector:temp_inspector_end)));
plot(imag(data(temp_inspector:temp_inspector_end)));
figure;
hold on;
plot(real(temp_data));
plot(imag(temp_data));
figure;
plot(abs(data(temp_inspector:temp_inspector_end)));
figure;
plot(unwrap(angle(data(temp_inspector:temp_inspector_end))));


% figure;
% pspectrum(data,sample_rate,'spectrogram');
st = 1 + 12.25 * 1024;
figure;
plot(unwrap(angle(data_tx(st:st + 1 * 1024))));
% plot(angle(data(start:start + 1.25 * 1024)));
% 
% figure;
% plot(unwrap(angle(down_chirp(1:1024)))); 

featureChirp = [up_chirp, up_chirp, down_chirp, down_chirp];
[corr, lag] = xcorr(simu_data, featureChirp);
% figure;
% plot(lag,abs(corr));

data_cutoff = data_tx(st:end);
y = fft(data_cutoff(1:1024) .* down_chirp);
figure;
plot(abs(y));
y = y(1:512);
pos = find(abs(y) == max(abs(y)));
bin = pos - 1;

f_start = BW / 2 - bin / 128 * BW ;
f_k = BW / chirp_duration;
t_flip = (f_start + BW / 2) / f_k;
chirp_conj = [chirp(0:1/sample_rate:t_flip - 1/sample_rate, f_start, chirp_duration, f_start - BW, 'linear',0,'complex') chirp(t_flip:1/sample_rate:chirp_duration - 1/sample_rate, f_start + BW, chirp_duration, f_start, 'linear',0,'complex')];

conj_sig = [down_chirp(1024 / 128 * bin + 1 : end) down_chirp(1: 1024/128*bin)];

st = 1 + 12.25 * 1024;

figure;
plot(unwrap(angle(chirp_conj))); 
figure;
plot(unwrap(angle(data_tx(st:st + 1024 - 1)))); 

figure;
plot(unwrap(angle(data_tx(st:st + 1024 - 1).* chirp_conj)));

figure;
plot(unwrap(angle(data_tx .* pkt_conj)));