clc
clear all
close all

set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);
linwidth = 2;

%%
% ===================================================================
%                         INPUT AUDIO SIGNAL
% ===================================================================

% 1. load audio file to MATLAB
[data, Fs] = audioread('B.wav'); % get sampling freq and data
playblocking(audioplayer(data, Fs));  
% sends audio signal y to the speaker at sample rate Fs
% Default Fs in MATLAB = 48000 Hz

% Signal in time Domain
L = length(data);      % Original Sample length
Ts = 1/Fs;             % Sampling time /s
t = (0:Ts:(L-1)/Fs);   % Time


% Plotting sound in time domain
figure(1)
subplot (2, 1, 1)
stem(t, data)
title('Audio Waveform in Time Domain - Before Filtering')
xlabel('Time/n')
ylabel('Amplitude')
grid on

% Using FFT to find signal in Frequency domain
 X = fft(data);          % original sample length
 n = pow2(nextpow2(L))   % Transform length 
 A1 = abs(X/n);          % Double array
 A2 = A1(1:L/2+1);       % L/2+1 when even and L/2+1.5 when odd
 A2(2:end-1) = 2*A2(2:end-1);
 f = ((0:1/n:(1-1/n))*Fs)  % frequency

% Plotting FFT Magnitude of Original Audio in Freq Domain
figure(1)
subplot(2,1,2)
plot(f(1:floor(n/2)), A1(1:floor(n/2))) %floor to nearest integer
title('Audio Waveform in Frequency Domain - Before FIltering')
xlabel('Frequency/Hz')
ylabel('Amplitude')
grid on

%%
% ===================================================================
%                   BUTTERWORTH LOWPASS FILTER
% ===================================================================

passband = 0.05*pi/pi; % convert to normalized frequency
stopband = 0.1*pi/pi;

passrip = -20*log10(0.9); % ripple in dB
stopatten = -20*log10(0.001); %stopband attenuation in dB

[Nb,Wnb] = buttord(passband,stopband,passrip,stopatten);
% lowest filter order = Nb
% cut-off freq vectors = Wnb
[Bb,Ab] = butter(Nb,Wnb);
% [b,a] = butter(n,Wn) returns the transfer function coefficients of an 
% nth-order lowpass digital Butterworth filter with normalized cutoff 
% frequency Wn

figure(2)
[Hb,W]= freqz(Bb,Ab,2048);
% [h,w] = freqz(b,a,n) returns the n-point frequency response vector h 
% and the corresponding angular frequency vector w for the digital filter 
% with transfer function coefficients stored in b and a.
subplot(3, 1, 1)
h = plot(W/pi,20*log10(abs(Hb)));
set(h,'LineWidth',linwidth);
title(['Frequeny Response: Butterworth Filter, Order = ',num2str(Nb)])
axis([0 1 -60, 5])
ylabel('Gain (dB)')
xlabel('Normalized Frequency (\times \pi rad/sample)')
grid on

figure(2)
[h1,t1]= impz(Bb,Ab);
% [h,t] = impz(b,a)  returns the impulse response of the digital filter 
% with numerator coefficients b and denominator coefficients a. 
% The function chooses the number of samples and returns the response 
% coefficients in h and the sample times in t.
subplot(3, 1, 2)
g = stem(t1, h1);
set(g,'LineWidth',linwidth);
title(['Impulse Response: Butterworth Filter, Order = ',num2str(Nb)])
ylabel('Amplitude')
xlabel('Samples / n')
grid on

figure(2)
[phi,w]= phasez(Bb,Ab,W);
% phasez(b,a,n) returns the n-point phase response vector phi and the 
% corresponding angular frequency vector w for the digital filter with the 
% transfer function coefficients stored in b and a.
subplot(3, 1, 3)
g = plot(phi, w);
set(g,'LineWidth',linwidth);
title(['Phase Response: Butterworth Filter, Order = ',num2str(Nb)])
ylabel('Phase Shift (rad)')
xlabel('Normalized Frequency (\times \pi rad/sample)')
grid on

figure(3)
j=zplane(Bb,Ab)     % Finding poles and zeroes
set(j,'LineWidth',linwidth);
title(['Pole-Zero Diagram: Butterworth Filter, Order = ',num2str(Nb)])
ylabel('Imaginary')
xlabel('Real')
grid on


%%
% ===================================================================
%                        FILTERING AUDIO SIGNAL
% ===================================================================

filtered_sound = filter(Bb, Ab, data);        %filter the audio

% Filtered audio in time domain
t1 = (0:Ts:(length(filtered_sound)-1)/Fs);
figure(4)
subplot(2, 1, 1)
stem(t1,filtered_sound)
title('Audio Waveform in Time Domain - After FIltering')
xlabel('Time/n')
ylabel('Amplitude')
grid on

% Filtered audio in freq domain
L1 = length(data);          % original sample length
n1 = pow2(nextpow2(L1));    % Transform length
X1 = fft(filtered_sound, n1);   % fft filter
A3 = abs(X1)/n1;
A4 = A3(abs(1:L/2+1));
A4(2:end-1)=2*A4(2:end-1);
f = ((0:1/n1:1-1/n1)*Fs);

% Plot filtered sound in freq domain
figure(4)
subplot(2, 1, 2)
plot(f(1:floor(n1/2)), A3(1:floor(n1/2))) %floor to nearest integer
title('Audio Waveform in Frequency Domain - After FIltering')
xlabel('Frequency/Hz')
ylabel('Amplitude')
grid on

%%
% ===================================================================
%                        OUTPUT AUDIO SIGNAL
% ===================================================================
playblocking(audioplayer(filtered_sound, 48000));
filename = 'B_filtered.wav';
audiowrite(filename,filtered_sound,48000);

