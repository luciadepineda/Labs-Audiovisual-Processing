function FreqResponse(h,Fs,Gain,color)
if nargin ==1
    Fs=1;
    Gain=1;
    color='b';
elseif nargin ==2
    Gain=1;
    color='b';
elseif nargin ==3    
    color='b';
end

[HI,WI]=freqz(h,1,512);
plot(WI/(2*pi)*Fs,20*log10(Gain*abs(HI)),color);
xlabel('Hz or Normalized frequency'); ylabel('Magnitude (dB)');
