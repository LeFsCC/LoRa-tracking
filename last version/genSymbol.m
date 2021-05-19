function sig = genSymbol(BW,SF,chirp_t,chirp_time,chirp_n,k,down)
    if down == 0
        f0 = -BW/2; % start freq
        f1 = BW/2;  % end freq
    else
        f0 = BW/2;
        f1 = -BW/2;
    end
    chirpI = chirp(chirp_t, f0, chirp_time, f1, 'linear', 90);
    chirpQ = chirp(chirp_t, f0, chirp_time, f1, 'linear', 0);
    baseline = complex(chirpI, chirpQ);
    baseline = repmat(baseline,1,2);
    clear chirpI chirpQ
    offset = round((2^SF - k) / 2^SF * chirp_n)
    sig = baseline(offset+(1:chirp_n));
end