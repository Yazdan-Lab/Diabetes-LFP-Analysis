Fs = 1250;

raw_pyr = LFP(:,chans(j+1));
pyr = BPfilter( raw_pyr,Fs,150,250);
pyr2 = pyr.^2; % filter for SWR (150-250 Hz) and square
s_pyr = smoothvect(pyr2, kernel);
Theta = BPfilter( raw_pyr,Fs,4,7);
Delta = BPfilter( raw_pyr,Fs,0.1,3);

H_Theta = abs(hilbert(Theta));
H_Delta = abs(hilbert(Delta));

S_Theta = smoothvect(H_Theta, kernel);
S_Delta = smoothvect(H_Delta, kernel);

Theta_Delta = S_Theta./S_Delta;
S_Theta_Delta = smoothvect(Theta_Delta, kernel);
%%
figure
for i = 1:625:length(raw_pyr)
    subplot(5,1,1)
    plot(raw_pyr(i:i+(1250*10)))
    
    subplot(5,1,2)
    plot(s_pyr(i:i+(1250*10)))
    
    subplot(5,1,3)
    plot(Theta(i:i+(1250*10)))
    
    subplot(5,1,4)
    plot(Delta(i:i+(1250*10)))
    
    subplot(5,1,5)
    plot(S_Theta_Delta(i:i+(1250*10)))
    hline(0.5)
    pause(0.5)
end