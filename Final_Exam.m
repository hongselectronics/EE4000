% Jeremy Hong
% 4/22/2015
% Wright State University
%
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
%
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
%
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
figure(2)
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
%
%%
% Problem 4: [15 points] Apply your FIR filter in the frequency domain and
% plot the magnitude of your results as 20*log10(abs(…)). Make sure that 
% ALL axes have proper values, limits and labels.
%
%%
% Problem 5: [20 points] Create an IIR lowpass filter with a cutoff of 
% 20 KHz (transition frequency 25 KHz) and keeps only the 10 KHz signal. 
% Make sure that the first filter sidelobe is at least 40 dB down. Plot 
% the normalized filter magnitude as 20*log10(abs(…)) and on a separate 
% plot show the unwrapped phase in degrees. Make sure that ALL axes have 
% proper values, limits and labels. Remember that the unwrap command has 
% input and output in radians.
%
%%
% Problem 6: [15 points] Apply your IIR filter in the frequency domain and 
% plot the magnitude of your results as 20*log10(abs(…)). Make sure that
% ALL axes have proper values, limits and labels.
%
%%
% Problem 7: [10 points] On a separate figure plot the results of applying 
% your FIR and IIR filters on the same axes in the frequency domain. Use 
% black for the FIR and red for the IIR. Plot the magnitude as 
% 20*log10(abs(…)) and on a separate figure show the phase. 
% Make sure that ALL axes have proper values, limits and labels.
%


