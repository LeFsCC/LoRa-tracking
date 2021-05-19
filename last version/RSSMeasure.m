
fid = fopen("../Data/4-2/back-forward-2",'rb');
raw_data =fread(fid,[2 inf],'float');
fclose(fid);
signal = complex(raw_data(1,:),raw_data(2,:));


mean_window_size = 100;
st = 1;
rssSeq = [];

while st + mean_window_size < length(signal)
    rssSeq = [rssSeq mean(abs(signal(st:st+mean_window_size)).^2)];
    st = st + mean_window_size;
end
rssSeq = [rssSeq mean(abs(signal(st:end)).^2)];
rssSeq = 10 * log2(rssSeq);
rssSeq = rssSeq(rssSeq(:) > -16);

width = 1000;

for i = 1:length(rssSeq)-width
    rssSeq(i) = mean(rssSeq(i:i+width));
end

figure;
plot(rssSeq, 'r', 'LineWidth', 3);
xlabel('time');
ylabel('dBm');


