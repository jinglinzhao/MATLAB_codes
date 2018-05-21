function [power,Y] = FUNCTION_iFFT(X)

n = length(X);
Y = ifft(X, n);
Y = fftshift(Y);    
power = abs(Y)*n; % squared? 