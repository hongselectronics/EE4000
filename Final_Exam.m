clc
clear
close all
% Jeremy Hong
% 4/22/2015
% Wright State University
% 
% 
% ©2014-2015 Hirsch Chizever (Linear Systems I/II) 
% ©2015 Hong's Electronics/Jeremy K Hong
%                               _       
%   /\  /\  ___   _ __    __ _ ( ) ___  
%  / /_/ / / _ \ | '_ \  / _` ||/ / __| 
% / __  / | (_) || | | || (_| |   \__ \ 
% \/ /_/   \___/ |_| |_| \__, |   |___/ 
%                        |___/     
%  _____  _              _                        _            
% |  ___|| |            | |                      (_)           
% | |__  | |  ___   ___ | |_  _ __   ___   _ __   _   ___  ___ 
% |  __| | | / _ \ / __|| __|| '__| / _ \ | '_ \ | | / __|/ __|
% | |___ | ||  __/| (__ | |_ | |   | (_) || | | || || (__ \__ \
% \____/ |_| \___| \___| \__||_|    \___/ |_| |_||_| \___||___/
%                                                             
% ©2015 Sierra Nevada Corporation
%    _____   _   _    _____ 
%   / ____| | \ | |  / ____|
%  | (___   |  \| | | |     
%   \___ \  | . ` | | |     
%   ____) | | |\  | | |____ 
%  |_____/  |_| \_|  \_____|
%
% Copyright notice & terms of use available at:
% https://github.com/hongselectronics/EE4000/blob/master/LICENSE.md
% View code revision history here:
% https://github.com/hongselectronics/EE4000/blob/master/Final_Exam.m
%
% Linear Systems II Final Exam: MATLAB Take Home Portion
% EE4000-01
% Spring 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Main problem definition: A signal generator puts out two sinusoids of 
% equal amplitude, one at 10 KHz and the other at 50 KHz. The signal is 
% sampled with an A/D operating at 200 KHz and collects data for 1 ms.
0;
% Definition and Precalculations:
% We have an Analog Devices A/D converter that operates at 200 KHz and 
% runs for 1 ms
fs = 200000;
%Calculate sample frequency
ts = 1/fs;
% time
t = 0:ts:.001;
% Define frequencies and sinusoids
% Our Tektronix signal generator outputs a signal at 10 KHz
Tek = 10000;
s1 = sin(2*pi*Tek*t);
% Our HP/Agilent/Keysight signal generator outputs a signal at 50 KHz
Key = 50000;
s2 = sin(2*pi*Key*t);
%%
% Problem 1 [10 points] Create a total signal in the time domain from the 
% two sinusoids and plot the real-valued signal. Make sure that ALL axes 
% have proper values, limits and labels.
0;
% total signal in time domain
s_t = s1 + s2;
% Plot signal in time domain
figure(1)
plot(t,s_t)
grid on
axis([0 .001 -2.5 2.5])
xlabel('Sample Time (sec)','fontweight','bold','fontsize',10)
% Since this is a signal generator we can assume that the y label is volts
ylabel('Amplitude (Volts)','fontweight','bold','fontsize',10)
title('10 KHz + 50 KHz Sinusoids Sampled @ 200 KHz','fontweight','bold','fontsize',10)

% 
%%
% Problem 2: [10 points] Compute the frequency domain of the total signal 
% using an 8192 FFT and plot your normalized magnitude as 20*log10(abs(…)) 
% Make sure that ALL axes have proper values, limits and labels.
%
nifft=8192;
fd=fftshift(fft(s_t,nifft));
fd=fd/max(fd);
xaxis=-fs/2:fs/nifft:fs/2-fs/nifft;
figure, grid, hold on
plot(xaxis,20*log10(abs(fd)),'k','linewidth',2)
xlabel('Frequency (Hz)','fontweight','bold','fontsize',10)
ylabel('Magnitude (dB)','fontweight','bold','fontsize',10)
title('Total Signal FFT','fontweight','bold','fontsize',10)
axis([-100000 100000 -100 0])

%%
% Problem 3: [20 points] Create an FIR lowpass filter with a cutoff of 
% 20 KHz and keeps only the 10 KHz signal. Make sure that the first filter 
% sidelobe is at least 40 dB down. Plot the normalized filter magnitude as 
% 20*log10(abs(…)) and on a separate plot show the unwrapped phase in
% degrees. Make sure that ALL axes have proper values, limits and labels.
% Remember that the unwrap command has input and output in radians.
0;
% Template from Hirsch Chizever's "CreateFIR_LowPass.m" example.
ntaps=301;
NormalizedFrequency=.20; % 1.0 equals fs/2
%
ZeroFreq=(ntaps+1)/2;
nfilter=NormalizedFrequency*(ntaps-1)+1; % filter indexes within passband
index=1:nfilter;
x=zeros(ntaps,1);
x(index)=1;
x=circshift(x,[ZeroFreq-round(nfilter/2) 0]); % Shift zero freq to center
x=x.*exp(-1i*2*pi*(0:ntaps-1)'*(ZeroFreq-1)/ntaps); % % Apply linear phase
y=ifft(x); % compute h(n)
ynew=real(y); 
xnew=fft(ynew,8192); % compute realized filter
%
% Plot realized log response versus windows
%
figure, grid, hold on
plot(linspace(-1,1-2/8192,8192),20*log10(abs(xnew)),'k','linewidth',2)
set(gca,'xtick',-1:.20:1)
set(gca,'fontweight','bold','fontsize',10)
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('Magnitude (dB)','fontweight','bold','fontsize',10)
title(['Low Pass Filter Response Cutoff=' num2str(NormalizedFrequency)],'fontweight','bold','fontsize',10)
set(gca,'fontweight','bold','fontsize',10)
axis([-1 1 -150 0])
%
% Plot realized phase response
%
figure, grid, hold on
plot(linspace(-1,1-2/8192,8192),unwrap(angle(xnew)).*180./(pi*100),'k','linewidth',2)
set(gca,'xtick',-1:.20:1)
set(gca,'fontweight','bold','fontsize',10)
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('Unwrapped Phase (deg)','fontweight','bold','fontsize',10)
title(['Low Pass Filter Response Cutoff=' num2str(NormalizedFrequency)],'fontweight','bold','fontsize',10)
set(gca,'fontweight','bold','fontsize',10)
axis([-1 1 -12000/100 5])
%%
% Problem 4: [15 points] Apply your FIR filter in the frequency domain and
% plot the magnitude of your results as 20*log10(abs(…)). Make sure that 
% ALL axes have proper values, limits and labels.
0;
% Template from Hirsch Chizever's "FilterTester.m" example.
N=301;
nifft=8192;
signal=sin(2*pi*(0:N-1)*.5)+sin(2*pi*(0:N-1)*.1);

taps=1000;
WinFunc='kaiser';
WinParm=8;
FilterType='lowpass';
band= 0.2;
% Using Hirsch Chizever's MakeFIR Function ©2014-2015 
[hn, Hw, w] = MakeFIR(FilterType, taps, band, nifft, 1, WinFunc,WinParm);
% 
filterSignal = ApplyFIR(Hw,signal,taps,nifft);
%
% Plot Apply FIR filter dB v. Normalized Frequency
%
xaxis=-1:2/nifft:1-2/nifft;
tdSignal=fftshift(ifft(signal, nifft));
tdSignal=tdSignal/max(tdSignal);
tdFilteredSignal=fftshift(ifft(filterSignal, nifft));
tdFilteredSignal=tdFilteredSignal/max(tdFilteredSignal);

figure, grid, hold on
plot(xaxis,20*log10(abs(tdSignal)),'k','linewidth',2)
plot(w,20*log10(abs(Hw)),'r','linewidth',2)
plot(xaxis,20*log10(abs(tdFilteredSignal)),'b','linewidth',2)
legend('Original Signal','Filter Response','Filtered Signal','location','south')
title(['Apply Low Pass FIR Filter Cutoff=' num2str(NormalizedFrequency)],'fontweight','bold','fontsize',10)
set(gca,'fontweight','bold','fontsize',10)
axis([-1 1 -250 0])
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('Magnitude (dB)','fontweight','bold','fontsize',10)
%%
% Problem 5: [20 points] Create an IIR lowpass filter with a cutoff of 
% 20 KHz (transition frequency 25 KHz) and keeps only the 10 KHz signal. 
% Make sure that the first filter sidelobe is at least 40 dB down. Plot 
% the normalized filter magnitude as 20*log10(abs(…)) and on a separate 
% plot show the unwrapped phase in degrees. Make sure that ALL axes have 
% proper values, limits and labels. Remember that the unwrap command has 
% input and output in radians.
0;
% Use Hirsch Chizever's MakeIIR Low Pass Filter template
%
% Make IIR lowpass filter

fc=.2; % Cutoff frequency
ft=.25; % Stop band

gamma1=.25;
gamma2=.75;

p1=gamma1; %DC Passband
p2=gamma2*exp(1i*pi*fc);
p3=gamma2*exp(-1i*pi*fc);

z1=exp(1i*pi*ft); % Stopband
z2=exp(-1i*pi*ft);
z3=-1; % fs/2

omega=-pi:2*pi/8192:pi-2*pi/8192; % Define omega
z=exp(1i*omega);

% Compute H(omega)
Hz=(z-z1).*(z-z2).*(z-z3)./(z-p1)./(z-p2)./(z-p3);
Klog=-16.0591;
Hz=Hz*10^(Klog/20);

% Plot magnitude of H(omega)
figure
plot(omega/pi,20*log10(abs(Hz)),'k','linewidth',2)
grid
axis([-1 1 -80 20])
set(gca,'fontweight','bold','fontsize',10)
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('|H(\Omega)|','fontweight','bold','fontsize',10)
title(['3rd Order Low Pass IIR Filter'],'fontweight','bold','fontsize',10)


% Plot phase of H(omega)
figure
plot(omega/pi,angle(Hz)*180/pi,'k','linewidth',2)
grid
axis([-1 1 -180 180])
set(gca,'fontweight','bold','fontsize',10)
xlabel('Normalized Frequency (\Omega/\pi)','fontweight','bold','fontsize',10)
ylabel('Normalized Angle (degrees)','fontweight','bold','fontsize',10)
title(['3rd Order Low Pass IIR Filter Phase Response'],'fontweight','bold','fontsize',10)

%%
% Problem 6: [15 points] Apply your IIR filter in the frequency domain and 
% plot the magnitude of your results as 20*log10(abs(…)). Make sure that
% ALL axes have proper values, limits and labels.
%
0;
% Plot magnitude of H(omega)
figure
hold on
plot(omega/pi,20*log10(abs(tdSignal)),'k','linewidth',2)
plot(omega/pi,20*log10(abs(tdSignal.*Hz)),'b','linewidth',2)
plot(omega/pi,20*log10(abs(Hz)),'r','linewidth',2)
grid
axis([-1 1 -150 0])
legend('Original Signal','Filter Response','Filtered Signal','location','south')
title(['Apply Low Pass IIR Filter Cutoff=' num2str(NormalizedFrequency)],'fontweight','bold','fontsize',10)
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('|H(\Omega)|','fontweight','bold','fontsize',10)
%%
% Problem 7: [10 points] On a separate figure plot the results of applying 
% your FIR and IIR filters on the same axes in the frequency domain. Use 
% black for the FIR and red for the IIR. Plot the magnitude as 
% 20*log10(abs(…)) and on a separate figure show the phase. 
% Make sure that ALL axes have proper values, limits and labels.
%
figure
hold on
grid on
plot(xaxis,20*log10(abs(tdFilteredSignal)),'k','linewidth',2)
plot(omega/pi,20*log10(abs(tdSignal.*Hz)),'r','linewidth',2)
axis([-1 1 -250 0])
legend('FIR Filtered Signal','IIR Filtered Signal','location','south')
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('Magnitude (dB)','fontweight','bold','fontsize',10)
title(['Low Pass FIR & IIR Filter Cutoff=' num2str(NormalizedFrequency)],'fontweight','bold','fontsize',10)

figure
hold on
grid on
plot(linspace(-1,1-2/8192,8192),unwrap(angle(xnew)).*180./(pi*100),'k','linewidth',2)
plot(omega/pi,angle(Hz)*180/pi,'r','linewidth',2)
xlabel('Normalized Frequency','fontweight','bold','fontsize',10)
ylabel('Normalized Angle (degrees)','fontweight','bold','fontsize',10)
legend('FIR Phase Response','IIR Phase Response','location','northeast')
title('Phase Responses of FIR & IIR Filters','fontweight','bold','fontsize',10)
