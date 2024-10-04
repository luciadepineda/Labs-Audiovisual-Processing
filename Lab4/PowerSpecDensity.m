function PowerSpecDensity(Signal,Fs,color)
if nargin==1
    Fs=1;
    color='b';
end
[Pxx,Wo]=periodogram(Signal,[ ],512);
Pxx=Pxx*2*pi/Fs; % Revisar reescalat Canvi variable eix W -> f
plot(Wo/(2*pi)*Fs,10*log10(Pxx),color);
grid on;
xlabel('Hz or Normalized frequency');
ylabel('PSD dB/Hz')
