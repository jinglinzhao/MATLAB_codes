f0 = 4;   %frequency of the sine wave
Fs = 140; %sampling rate
Ts = 1/Fs; %sampling time interval
t = 0:Ts:2.4-Ts; %sampling period
n = length(t); %number of samples
y = 2*sin(2*pi*f0*t) ; %+ 2*sin(2*pi*f1*t) ; %the sine curve
 
%plot the cosine curve in the time domain
figure;
subplot(5,1,1);
%sinePlot = figure;
plot(t,y)
xlabel('time (seconds)')
ylabel('y(t)')
title('Sample Sine Wave')
grid



[YfreqDomain_half, YfreqDomain, frequency_half, frequency] = positiveFFT(y,Fs);



%figure;
subplot(5,1,2);
stem(frequency_half,abs(YfreqDomain_half));
xlabel('Freq (Hz)')
ylabel('Amplitude')
title('Using the positiveFFT function')
grid on
% axis([0,100,0,3])

subplot(5,1,3);
stem(frequency,abs(YfreqDomain));