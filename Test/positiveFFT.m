function [Y_freq_half,Y_freq,freq_half,freq]=positiveFFT(x,Fs)
   N=length(x);   %get the number of points
   k=0:N-1;        %create a vector from 0 to N-1
   T=N/Fs;          %get the frequency interval
   freq=k/T;     %create the frequency range
   %X=fft(x)/N*2;  % normalize the data
   Y_freq=fft(x)/N;
 
   %only want the first half of the FFT, since it is redundant
   cutOff = ceil(N/2);  
   %take only the first half of the spectrum
    freq_half = freq(1:cutOff);
   Y_freq_half= 2*Y_freq(1:cutOff);