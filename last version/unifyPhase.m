function [ sig2 ] = unifyPhase( signal,Fs,fsquare, chirp_n, payload, offset)
% 先变成正弦波，再让相位连续
% 调整每个chirp起始位置的相位连续性

sig1=signal(1:chirp_n);
for k=1:length(signal)/chirp_n-1
    b1 = signal(k*chirp_n+1:(k+1)*chirp_n);
    sig1=[sig1,b1];
    new_b1=axisConversion(b1,(angle(sig1(k*chirp_n+1))-angle(sig1(k*chirp_n)))-2*pi/(Fs/fsquare));
    ampli=abs(b1);
    anglNew=angle(b1)-((angle(sig1(k*chirp_n+1))-angle(sig1(k*chirp_n)))-2*pi/(Fs/fsquare));
    for n=1:length(b1)
    new_b1(n)=ampli(n)*(cos(anglNew(n))+i*sin(anglNew(n)));
    end
    
    sig1(k*chirp_n+1:(k+1)*chirp_n)=new_b1;
end

if payload == 1
% 调整每个chirp中间offset导致的相位连续性
    offset;
    sig2=sig1(1:offset);
    for m=1:length(signal)/chirp_n-1
        c1=sig1(offset+chirp_n*(m-1)+1:offset+chirp_n*m);   
        sig2=[sig2,c1];
        new_c1=axisConversion(c1,(angle(sig2(offset+chirp_n*(m-1)+1))-angle(sig2(offset+chirp_n*(m-1))))-2*pi/(Fs/fsquare));
        sig2(offset+chirp_n*(m-1)+1:offset+chirp_n*m)=new_c1;
    end
    % 调整最后一段的相位
    d1=sig1(offset+chirp_n*(length(signal)/chirp_n-1)+1:end);
    new_d1=axisConversion(d1,(angle(d1(1))-angle(sig2(end)))-2*pi/(Fs/fsquare));
    sig2=[sig2,new_d1];
end

if payload == 0
    sig2=sig1;
end

end

