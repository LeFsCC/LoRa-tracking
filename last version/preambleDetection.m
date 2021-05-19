function [ position ] = preambleDetection( x,upChirp,downChirp,chirp_n )

featureChirp =[upChirp,downChirp];
[corr, lag] = xcorr(x, featureChirp);
figure
plot(abs(corr))
cLag=find(abs(corr)==max(abs(corr)))
position=lag(cLag);

%corrThresh = max(abs(corr))/1.1;
%cLag = find(abs(corr) > corrThresh);
%select_cLag=[];
%for k=1:length(cLag)
    %if find(cLag==cLag(k)+chirp_n)&find(cLag==cLag(k)+2*chirp_n)&find(cLag==cLag(k)+3*chirp_n)&find(cLag==cLag(k)+4*chirp_n)&find(cLag==cLag(k)+5*chirp_n)&find(cLag==cLag(k)+6*chirp_n)&find(cLag==cLag(k)+7*chirp_n)
    %if find(cLag==cLag(k)+chirp_n)&find(cLag==cLag(k)+2*chirp_n)&find(cLag==cLag(k)+3*chirp_n)
        %select_cLag=[select_cLag,cLag(k)];
    %end
%end

%Index= select_cLag(find(abs(corr(select_cLag))==max(abs(corr(select_cLag))))); %
%position=abs(lag(Index));


end




