n=[0:15999];
    Noise1=randn(size(n));
    SP1(i, :) = Noise1;
Fs=8000;
figure
plot(n/Fs,Noise1);
title('Process 1')

%% 

SP1=[];
i=1;
while i <= 2000
    n=[0:15999];
    Noise1=randn(size(n));
    SP1(i, :) = Noise1;
    i=i+1;
end

figure
subplot(4,1,1);
plot(SP1(1,:));
title('Realization 1')
subplot(4,1,2);
plot(SP1(2,:));
title('Realization 2')
subplot(4,1,3);
plot(SP1(3,:));
title('Realization 3')
subplot(4,1,4);
plot(SP1(4,:));
title('Realization 4')
sgtitle('Realizations stochastic process 1')

%% 
ni=983
figure
hist(SP1(:,ni),100);
DeltaX=(max(SP1(:,ni))-min(SP1(:,ni)))/100;
K1=DeltaX*size(SP1,1);
GaussPDF=normpdf([-4:0.01:4],0,1); % Exact pdf
K2=0.01*sum(GaussPDF);
hold on
plot([-4:0.01:4],GaussPDF*K1/K2,'r')
title('Histogram instant 983')


%% 
n2=[0:15999];
A=exp(-(n2/8000-1).^2/0.1);
Noise2=sqrt(A).*randn(size(A));
Fs=8000;
figure
plot(n/Fs,Noise2);
title('Process 2');

%% 

SP2=[];
i=1;
while i <= 2000
    n2=[0:15999];
    A=exp(-(n2/8000-1).^2/0.1);
    Noise2=sqrt(A).*randn(size(A));
    SP2(i, :) = Noise2;
    i=i+1;
end

%% 
figure
subplot(4,1,1);
plot(SP2(1,:));
title('Realization 1')
subplot(4,1,2);
plot(SP2(2,:));
title('Realization 2')
subplot(4,1,3);
plot(SP2(3,:));
title('Realization 3')
subplot(4,1,4);
plot(SP2(4,:));
title('Realization 4')
sgtitle('Realizations stochastic process 2')
%% 
figure
hist(SP2(:,5000),100);
DeltaX=(max(SP2(:,5000))-min(SP2(:,5000)))/100;
K1=DeltaX*size(SP2,1);
GaussPDF=normpdf([-4:0.01:4],0,1); % Exact pdf
K2=0.01*sum(GaussPDF);
hold on
plot([-4:0.01:4],GaussPDF*K1/K2,'r')
title('Histogram sample 5000')

%% 
h1=ones(1,10);
h1=h1./(sqrt(h1*h1'));
n3=[0:15999];
NoiseX=randn(size(n3));
NoiseY=conv(NoiseX,h1,'same');
figure
plot(n3/Fs,NoiseY);
title('Process 3')

%% 

SP3=[];
i=1;
while i <= 2000
    h1=ones(1,10);
    h1=h1./(sqrt(h1*h1'));
    n3=[0:15999];
    NoiseX=randn(size(n3));
    NoiseY=conv(NoiseX,h1,'same');
    SP3(i, :) = NoiseY;
    i=i+1;
end

%% 
[Rx1 Tau1]=IPALab2_XCorr(SP1,127)
plot(Tau1,Rx1)
title('Autocorrelation SP1')
%% 
fft1=abs(fft(Rx1));
plot(fft1);
ylim([0 10])
title('Power spectral density SP1');
%% 
fft3=abs(fft(Rx3));
plot(fft3);
ylim([0 10])
title('Power spectral density SP3');
%% 
SP4=[];
i=1;
while i <= 2000
    h1=ones(1,100);
    h1=h1./(sqrt(h1*h1'));
    n3=[0:15999];
    NoiseX=randn(size(n3));
    NoiseY=conv(NoiseX,h1,'same');
    SP4(i, :) = NoiseY;
    i=i+1;
end

figure
subplot(4,1,1);
plot(SP4(1,:));
title('Realization 1')
subplot(4,1,2);
plot(SP4(2,:));
title('Realization 2')
subplot(4,1,3);
plot(SP4(3,:));
title('Realization 3')
subplot(4,1,4);
plot(SP4(4,:));
title('Realization 4')
sgtitle('Realizations stochastic process 4')
%% 
h1=ones(1,100);
h1=h1./(sqrt(h1*h1'));
n3=[0:15999];
NoiseX=randn(size(n3));
NoiseY=conv(NoiseX,h1,'same');
Fs=8000;
figure
plot(n3/Fs,NoiseY);
title('Process 4')

%% 
[Rx4 Tau4]=IPALab2_XCorr(SP4,127);
plot(Tau4,Rx4);
title('Autocorrelation SP4');
%% 
fft4=abs(fft(Rx4));
plot(fft4);
ylim([0 10]);
title('Power spectral density SP4');
%% 
m2=mean(SP2(728,:))
v2=var(SP2(728,:))

%% 
figure
subplot(3,1,1);
plot(Power(2,:));
title('Realization 2')
subplot(3,1,2);
plot(Power(4,:));
title('Realization 4')
subplot(3,1,3);
plot(Power(7,:));
title('Realization 7')
sgtitle('Realizations stochastic process')

%% 

[y Fs]=Audio_Read('Audio1.wav');

%% 

[Rxy Tauy] = IPALab2_XCorr(y, 2047)
plot(Tauy,Rxy);
title('Autocorrelation of the process');

%% 

ffty = abs(fft(Rxy));
plot(ffty);
title('Power spectral density of the process');


