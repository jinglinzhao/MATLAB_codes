function [fshift,power,Y] = FUNCTION_FFT(X, Fs)

% To use the fft function to convert the signal to the frequency domain, 
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.

% n       = 2^nextpow2(length(X));
n = length(X);
fshift  = Fs*((0.5:n-0.5)/n-0.5);
% fshift  = transpose(fshift);

% Power is the squared magnitude of a signal's Fourier transform,
% normalized by the number of frequency samples.
Y = fft(X, n);
Y = fftshift(Y);    
power = abs(Y).^2/n; % squared? 