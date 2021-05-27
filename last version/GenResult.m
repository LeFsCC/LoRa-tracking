BW = 500000;
SF = 11;
sample_rate = 1e6;
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
pkt_size = 12.25 + 13;
move_rate = 3e8 / 902e6;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');
file_name = "2-3";

fid = fopen(strcat("../Data/1D Back forward/", file_name),'rb');
raw_data =fread(fid,[2 inf],'float');

fclose(fid);
raw_data_complex = complex(raw_data(1,:),raw_data(2,:));
% raw_data_complex = raw_data_complex(1:0.93*end);

[data, flip_idx, position, sin_signal,index_ls] = PacketContentDetect(raw_data_complex,BW,SF,sample_rate,pkt_size);
[cfc, diff2, diff1] = CurveFitRecover(sin_signal,flip_idx,index_ls,chirp_samples, sample_rate, BW);

figure;
plot(unwrap(angle(diff1)));
figure;
plot(unwrap(angle(diff2)));

f1 = strcat('../Data/5-22-result/', file_name, '-diff1.mat');
f2 = strcat('../Data/5-22-result/', file_name, '-diff2.mat');
save(f1, 'diff1');
save(f2, 'diff2');



