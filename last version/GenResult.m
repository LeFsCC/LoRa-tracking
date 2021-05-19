BW = 500000;

SF = 11;
sample_rate = 1e6;
chirp_duration = 2^SF / BW;
chirp_samples = 2^SF * sample_rate / BW ;
pkt_size = 12.25 + 13;
move_rate = 3e8 / 902e6;
down_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, BW/2, chirp_duration, -BW/2,'linear',0,'complex');
up_chirp = chirp(0:1/sample_rate:chirp_duration - 1/sample_rate, -BW/2, chirp_duration, BW/2,'linear',0,'complex');
fid = fopen("../Data/3-30/interfare-2",'rb');
raw_data =fread(fid,[2 inf],'float');

fclose(fid);
raw_data_complex = complex(raw_data(1,:),raw_data(2,:));

[data, flip_idx, position, sin_signal,index_ls] = PacketContentDetect(raw_data_complex,BW,SF,sample_rate,pkt_size);

[cfc, diff2, diff1] = CurveFitRecover(sin_signal,flip_idx,index_ls,chirp_samples, sample_rate, BW);

figure;
set(0,'defaultfigurecolor','w');
plot(unwrap(angle((diff1))));
title('目标物体反射信号相位');
xlabel('chirp number');
ylabel('phase/rad');
average_move = BreathFilter(diff1);
average_move = average_move(1:100:length(average_move));
intersection_tracking(average_move);


