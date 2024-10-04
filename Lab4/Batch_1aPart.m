%% 1.1 To + Soroll. Rx coneguda. hopt =Rx^-1·rd

clear all
close all

% Parameters
N=10000;
Q=50;
F0=1/16;
D=1;
Sigma2=1;
P=0.5;

% Signal and interference generation
w=sqrt(Sigma2)*randn(N,1);
Sigma2=std(w).^2;  % We update the real value for Sigma2 (approx. 1) 

Omega=2*pi*F0;
x=sqrt(2*P)*cos(Omega*[0:N-1]'+2*pi*rand);

v=w+x;

% Calculate correlation and R and rd vectors

m=[0:Q]';
rx=zeros(length(m),1);
rx(1)=Sigma2;
rx=rx+P*cos(Omega*m);

Rx  = toeplitz(rx(1:end-1),rx(1:end-1));
rxd = rx(2:end);

% Filter design
hopt1=Rx\rxd;
aopt1=[1;zeros(D-1,1);-hopt1];

% Wiener Filter frequency response
figure(1)
FreqResponse(hopt1)
hold on 
PowerSpecDensity(v,1,'r')
axis([0 0.5 -50, 40])
legend('Wiener filter Freq Response','Power Spectral Density Input signal')

%% 

% Global Filter frequency response
figure(2)
FreqResponse(aopt1,1,1,'k')
hold on 
PowerSpecDensity(v,1,'r')
axis([0 0.5 -15, 40])
legend('Global filter Freq Response','Power Spectral Density Input signal')

%% 

% Filter signal 
xestimate=filter([zeros(D,1);hopt1],1,v);
westimate=v-xestimate;

%Plot signals
figure(3)
subplot(3,1,1)
plot(x(1:100),'r');hold on
plot(v(1:100),'b')
title('x[n] (red) and v[n] (blue)')
subplot(3,1,2)
plot(x(1:100),'r');hold on
plot(xestimate(1:100),'mo-')
title('x[n] (red) and x estimated (pink)')
subplot(3,1,3)
plot(w(1:100),'r');hold on
plot(westimate(1:100),'mo-')
title('w[n] (red) and w estimated (pink)')


%% 1.2.1 To + Soroll.  Resolució amb LMS

clear all
close all

rng(1)
% Parameters
N=10000;
Q=32;
F0=1/16;
D=1;
mu=0.004;
P=0.5;
Sigma2=1;

% Signal and interference generation
w=sqrt(Sigma2)*randn(N,1);
x=sqrt(2*P)*cos(2*pi*F0*[0:N-1]'+2*pi*rand);
v=w+x;


% LMS Algorithm
rx_0_Estimate=v'*v/length(v);
disp(sprintf('Reference to set the mu constant:%1.4f',2/(Q*rx_0_Estimate)))

mu=0.03
h=zeros(Q,length(v)-2);
e=zeros(1,length(v)-2);
y=zeros(1,length(v)-2);

h(:,1)=ones(Q,1);                  % Arbitrary Initialization
for nn=Q+1:length(v)-D
    xd=v(nn-1:-1:nn-Q);
    y(nn-Q)=h(:,nn-Q)'*xd;
    e(nn-Q)=v(D+nn-1)-y(nn-Q);     % e -> The error is the cleaned signal.
    h(:,nn-Q+1)=h(:,nn-Q)+mu*(xd*e(nn-Q));
end
xestimate=y;            
westimate=e;

% Convergence of the coefficients
figure(11)
plot(h(1,:),'b');hold on;plot(h(2,:),'g')
title('Convergence of the two first coeff. of the filter h[0] and h[1]')

% Wiener Filter frequency response
figure(13)
hcoef=h(:,3*N/4);
FreqResponse(hcoef)
hold on 
PowerSpecDensity(v,1,'r')
PowerSpecDensity(xestimate,1,'g')
axis([0 0.5 -70, 40])
legend('Wiener filter Freq Response','Power Spectral Density Input signal','Power Power Spectral Density Output Wiener Filter')

% Global Filter frequency response
figure(14)
acoefs=[1;zeros(D-1,1);-h(:,3*N/4)];
FreqResponse(acoefs,1,1,'k')
hold on 
PowerSpecDensity(v,1,'r')
axis([0 0.5 -70, 40])
legend('Global filter Freq Response','Power Spectral Density Input signal')


figure(15)
IIni=N/2; IEnd=N/2+100;
subplot(3,1,1)
 plot(x(IIni:IEnd),'r');hold on
 plot(v(IIni:IEnd),'b')
subplot(3,1,2)
 plot(x(IIni:IEnd),'r');hold on
 plot(xestimate(IIni-D+1:IEnd-D+1),'mo-')
subplot(3,1,3)
 plot(w(IIni:IEnd),'r');hold on
 plot(westimate(IIni-D+1:IEnd-D+1),'mo-')


%% 1.2.2 To + Veu. Resolució amb LMS

clear all
close all

% Parameters
Q=32;
F0=1/16;
D=1;
mu=0.01;
P=0.5;


% Signal and interference generation
[w,Fs]=Audio_Read('Force.mp3');   % Voice
N=length(w);
Nf=N-D-Q+1;

x=sqrt(2*P)*cos(2*pi*F0*[0:N-1]'+2*pi*rand);
% Interference Chirp
% x=sqrt(2*P)*chirp([0:N-1]'/Fs,0,(N-1)/Fs,0.5*Fs);

v=w+x;

% LMS Algorithm
rx_0_Estimate=v'*v/length(v);
disp(sprintf('Reference to set the mu constant:%1.4f',2/(Q*rx_0_Estimate)))

h=zeros(Q,length(v)-2);
e=zeros(1,length(v)-2);
y=zeros(1,length(v)-2);

h(:,1)=ones(Q,1);                  % Arbitrary Initialization
for nn=Q+1:length(v)-D
    xd=v(nn-1:-1:nn-Q);
    y(nn-Q)=h(:,nn-Q)'*xd;
    e(nn-Q)=v(D+nn-1)-y(nn-Q);     % e -> the error is the Voice signal cleaned.
    h(:,nn-Q+1)=h(:,nn-Q)+mu*(xd*e(nn-Q));

% Uncomment to change the tone frequency after N/2 iterations 
%   Emulates a non-stationary process
%   See how the filter tracks the new frequency.
    if nn==N/2
        Omega=2*pi*1/8;
        x(N/2:end)=sqrt(2*P)*cos(Omega*[N/2-1:N-1]'+2*pi*rand);
        v=w+x;
    end
%      
end
xestimate=y;
westimate=e;

% Convergence of the coefficients
figure(21)
plot(h(1,:),'b');hold on;plot(h(2,:),'g')
title('Convergence of the two first coeff. of the filter h[0] and h[1]')

% PSD dels senyals
figure(22)
PowerSpecDensity(v,1,'r')
hold on
PowerSpecDensity(w,1,'k')
PowerSpecDensity(x,1,'g')
axis([0 0.5 -100, 40])
legend('PSD Veu + Sinusoide','PSD Veu','PSD Sinusoide')


% Wiener Filter frequency response
% Avaluem el filtre en n=10000 instant del senyal de veu té molta potència
figure(23)
hcoef=h(:,10000);
FreqResponse(hcoef)
hold on 
PowerSpecDensity(v,1,'r')
PowerSpecDensity(xestimate,1,'g')
axis([0 0.5 -70, 40])
legend('Wiener filter Freq Response','Power Spectral Density Input signal','Power Power Spectral Density Output Wiener Filter')

% Global Filter frequency response
figure(24)
acoefs=[1;zeros(D-1,1);-h(:,10000)];
FreqResponse(acoefs,1,1,'k')
hold on 
PowerSpecDensity(v,1,'r')
axis([0 0.5 -70, 40])
legend('Global filter Freq Response','Power Spectral Density Input signal')


figure(25)
IIni=N/2-50; IEnd=N/2+250;subplot(3,1,1)
subplot(3,1,1)
 plot(x(IIni:IEnd),'r');hold on
 plot(v(IIni:IEnd),'b')
 title('x[n] (red) and v[n] (blue)')
subplot(3,1,2)
 plot(x(IIni:IEnd),'r');hold on
 plot(xestimate(IIni-D+1:IEnd-D+1),'mo-')
 title('x[n] and x estimated')


 
% Sound signals
% soundsc(w,Fs);
% soundsc(v,Fs);
% soundsc(xestimate(100:end),Fs); % We eliminate the transitory
soundsc(westimate(100:end),Fs); % We eliminate the transitory


% Plots the filter frequency response each 100 iterations to see how the filter converges 
%
 figure(33)
 PowerSpecDensity(v,1,'r')
% 
 for nn=1:size(h,2)/10000
    figure(33)
    hold on 
    FreqResponse(h(:,nn*10000))
    FreqResponse([1;zeros(D-1,1);-h(:,nn*10000)],1,1,'k')
    hold off
    axis([0 0.5 -70, 40])
    legend('Power Spectral Density Input signal','Wiener filter Freq Response','Global filter Freq Response')
     
    pause
    a=get(gca,'children');
    delete(a(1)) 
    delete(a(2))
 end
